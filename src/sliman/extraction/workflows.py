"""
sliman/extraction/workflows.py

Dylan Ross (dylan.ross@pnnl.gov)

    This module houses functions that perform high level data extraction workflows
"""

import json
from time import time
from tempfile import TemporaryDirectory
from typing import Optional, Tuple, List, Dict
import sqlite3
import os

import numpy as np
from matplotlib import pyplot as plt

from mzapy.peaks import find_peaks_1d_gauss, calc_gauss_psnr, find_peaks_1d_localmax, _gauss
from mzapy.view import plot_atd, plot_spectrum, add_peaks_to_plot

from sliman._util import tol_from_ppm
from sliman.extraction._uimf import UReader, UReaderSetMzParams
from sliman.extraction._mobie import MReader
from sliman.extraction._targeted import (
    select_survey_atd_sparse, 
    select_full_atd_sparse, 
    select_full_spectrum_sparse
)


# worker threads for spectrum decoding
_UREADER_WORKERS = 8

# size of window around target arrival time to extract if
# using a static window
STATIC_WINDOW_SCANS = 500

# multiplier for peak FWHM that determines how far data 
# is extracted on either side of a fitted peak
FWHM_MULT = 2


def _survey_analysis(rdr, 
                     target_labels, 
                     target_mzs, 
                     mz_ppm, 
                     frame_mod, 
                     scan_mod, 
                     i_min, 
                     scan_min, 
                     scan_max,
                     min_peak_ht, 
                     min_peak_fwhm, 
                     max_peak_fwhm, 
                     min_dist
                     ):
    """
    Extract survey data, try to find peaks for targets,
    return set of targets to include for full data extraction based on peaks that were found

    Parameters
    ----------
    rdr : ``sliman.extraction._uimf._UReader``
        interface for reading data from UIMF files
    target_labels : ``list(str)``
        list of labels for targets
    target_mzs : ``list(float)``
        list of target m/zs
     mz_ppm : ``float``
        m/z tolerance for ATD selection in ppm
    frame_mod : ``int``
        determines how densely the frames are sampled in survey data
        1 in {frame_mod} frames are selected
    scan_mod : ``int``
        determines how densely the arrival time (scans) is sampled in survey data
        1 in {scan_mod} scans are selected
    i_min : ``float``
        minimum intensity cutoff for individual spectrum points
    min_peak_ht : ``float``
    min_peak_fwhm : ``float``
    max_peak_fwhm : ``float``
        ATD peak fitting parameters
    min_dist : ``float``
        minimum distance (in arrival time) between peaks

    Returns
    -------
    found_labels : ``list(str)``
    found_mzs : ``list(float)``
    found_dts : ``list(float)``
    found_fwhms : ``list(float)``
        Lists of targets with associated ATD peak drift time and FWHM for peaks that were found
    """
    found_labels, found_mzs, found_dts, found_fwhms = [], [], [], []
    print('extracting survey data ... ', end='')
    t0 = time()
    survey_data = rdr.extract_survey_data_sparse(frame_mod, scan_mod, i_min, scan_min, scan_max)
    print('done ({:.1f} s)'.format(time() - t0))
    for lbl, mz in zip(target_labels, target_mzs):
        print(lbl, mz, end=' -> ')
        mz_tol = tol_from_ppm(mz, mz_ppm)
        s, i = select_survey_atd_sparse(mz, mz_tol, max(rdr.scans_per_frame), scan_mod, survey_data)
        # look for up to 2 peaks
        pkdts, pkhts, pkwts = find_peaks_1d_gauss(s, i, 0.1, min_peak_ht, min_peak_fwhm, max_peak_fwhm, 3, True)
        # filter peaks so multiple peaks with about the same arrival time are not kept
        unique_peaks = []
        for pkdt, pkht, pkwt in zip(pkdts, pkhts, pkwts):
            unq = True
            for upk in unique_peaks:
                if abs(pkdt - upk[0]) < min_dist:
                    unq = False
                    break
            if unq:
                unique_peaks.append((pkdt, pkht, pkwt))
        found = False
        for pkdt, pkht, pkwt in unique_peaks:
            found_labels.append(lbl)
            found_mzs.append(mz)
            found_dts.append(pkdt)
            found_fwhms.append(pkwt)
            found = True
            print('arrival time: {:.0f} +/- {:.0f} scans ({:.3e})'.format(pkdt, pkwt, pkht))
        if not found:
            print('no ATD peak found')
        else :
            # plot atd
            ax = plot_atd(s, i, figsize=(3, 2), dt_unit='scans',
                          dt_range=(min(pkdts) - max(pkwts) * FWHM_MULT, max(pkdts) + max(pkwts) * FWHM_MULT))
            for pkdt, pkht, pkwt in unique_peaks:
                ax.plot(pkdt, pkht, 'kx', ms=7, mew=1.5)
            plt.show()
            plt.close()
    return found_labels, found_mzs, found_dts, found_fwhms


"""
FULL ANALYSIS RESULTS 

results = {
    "R.QTLLLRPGGK.W_z2": [
        # ATD peak 1
        {
            "mz": 541.8379, 
            "pk_params": (
                12393.3,          # arrival time (scans)
                98255,            # peak height
                150.0,            # FWHM 
                17.3              # SNR
            ), 
            "atd": (
                np.array([...]),  # ATD arrival time (scans)
                np.array([...])   # ATD intensities
            ), 
            "spectrum": (
                np.array([...]),  # MS1 m/z (M-1.5 → M+2.5)
                np.array([...])   # MS1 intensities
            ),
        },
        # ATD peak 2 (if second ATD peak found)
        {...} 
    ],
    "K.GLVVDMDGFEEERK.L_z3|K.NPNTSEPQHLLVM*K.G_z3": [...],
    ...
}
"""


def _full_analysis(rdr, 
                   survey_targets, 
                   i_min, 
                   mz_ppm, 
                   min_peak_ht, 
                   min_peak_fwhm, 
                   max_peak_fwhm, 
                   full_min_psnr,
                   static_window
                   ):
    """
    """
    found_ = {}
    # determine the scans that need to be included in full data extraction
    include_scans = set()
    for dtf, wtf in zip(survey_targets[2], survey_targets[3]):
        if static_window: 
            for s in range(int(dtf - STATIC_WINDOW_SCANS), int(dtf + STATIC_WINDOW_SCANS) + 1):
                include_scans.add(s)
        else:
            for s in range(int(dtf - FWHM_MULT * wtf), int(dtf + FWHM_MULT * wtf)):
                include_scans.add(s)
    print('{} scans included in full data extraction'.format(len(include_scans)))
    print('extracting full data ... ', end='')
    t0 = time()
    full_data = rdr.extract_full_data_sparse(i_min, include_scans=include_scans)
    print('done ({:.1f} s)'.format(time() - t0))
    for lbl, mz in zip(survey_targets[0], survey_targets[1]):
        print(lbl, mz, end=' -> ')
        mz_tol = tol_from_ppm(mz, mz_ppm)
        atd = select_full_atd_sparse(mz, mz_tol, max(rdr.scans_per_frame), full_data)
        # look for up to 2 peaks
        pkdts, pkhts, pkwts = find_peaks_1d_gauss(*atd, 0.1, min_peak_ht, min_peak_fwhm, max_peak_fwhm, 2, True)
        found = False
        for pkdt, pkht, pkwt in zip(pkdts, pkhts, pkwts):
            psnr = calc_gauss_psnr(*atd, (pkdt, pkht, pkwt))
            if psnr >= full_min_psnr:
                found_data = {
                    'mz': mz,
                    'pk_params': (pkdt, pkht, pkwt, psnr),
                    'atd': atd,
                    'spectrum': select_full_spectrum_sparse(mz, 0.025, pkdt - pkwt, pkdt + pkwt, full_data),
                }
                if not found:
                    found_[lbl] = [found_data]
                else:
                    found_[lbl].append(found_data)
                found = True
                print('arrival time: {:.0f} +/- {:.0f} scans ({:.3e}, {:.1f})'.format(pkdt, pkwt, pkht, psnr))
        if not found:
            print('no ATD peak found')
        else:
            # plot atd
            pkdt1, pkwt1 = pkdts[0], pkwts[0]
            ax = plot_atd(*atd, figsize=(3, 2), dt_unit='scans',
                          dt_range=(pkdt1 - pkwt1 * FWHM_MULT, pkdt1 + pkwt1 * FWHM_MULT))
            add_peaks_to_plot(ax, pkdts, pkhts, pkwts, c='k')
            plt.show()
            plt.close()
            # plot MS1 spectrum
            ax = plot_spectrum(*found_[lbl][0]['spectrum'],
                               figsize=(3, 2))
            ax.axvline(mz, c='k', ls='--', lw=0.75)
            plt.show()
            plt.close()
    return found_
        

def _target_is_in_selected_dataset(sel_dslbls: List[str], 
                                   dsid_to_dslbl: Dict[int, str],
                                   target_dsids: str
                                   ) -> bool :
    """
    Determines whether the set of dataset ids for a target contains at least one of the datasets
    that were selected for analysis. Uses dsid_to_dslbl to map dataset ids to dataset labels.
    Returns a bool
    """
    for tlbl in [dsid_to_dslbl[int(dsid)] for dsid in target_dsids.split(",")]:
        if tlbl in sel_dslbls:
            return True
    return False


def extract(uimf: str, 
            amt_tag_db: str,
            sel_dsets: List[str],
            min_tmz: float,
            max_tmz: float,
            mz_ppm: float, 
            frame_mod: int, 
            scan_mod: int,
            min_peak_ht: float, 
            min_peak_fwhm: float, 
            max_peak_fwhm: float,
            min_dist: float, 
            i_min: float, 
            full_min_psnr: float,
            mz_cal_params: Optional[Tuple[float, float]] = None
            ) :
    """
    Parameters
    ----------
    uimf : ``str``
        path to UIMF file
    amt_tag_db : ``str``
        AMT tag database 
    sel_dsets : ``list(str)``
        list of selected dataset labels to get targets from
    min_tmz : ``float``
    max_tmz : ``float``
        min/max m/z range for targets
    mz_ppm : ``float``
        m/z tolerance for ATD selection, in ppm
    frame_mod : ``int``
        determines how densely the frames are sampled in survey data
        1 in {frame_mod} frames are selected
    scan_mod : ``int``
        determines how densely the arrival time (scans) is sampled in survey data
        1 in {scan_mod} scans are selected
    min_peak_ht : ``float``
    min_peak_fwhm : ``float``
    max_peak_fwhm : ``float``
        ATD peak fitting parameters
    min_dist : ``float``
        minimum distance (in arrival time) between peaks from survey analysis
    i_min : ``float``
        minimum intensity cutoff for individual spectrum points
    full_min_psnr : ``float``
        minimum peak SNR to accept peak from full data
    mz_cal_params : ``tuple(float, float)``, optional
        apply the m/z calibration using the parameters (slope, intercept) if provided

    Returns
    -------
    extraction results
    """
    assert os.path.isfile(amt_tag_db), f"AMT tag database {amt_tag_db} not found"
    amt_con = sqlite3.connect(amt_tag_db)
    amt_con.row_factory = sqlite3.Row  # enable accessing query result row elements by key
    amt_cur = amt_con.cursor()
    if mz_cal_params is not None:
        mz_cal_slope, mz_cal_intercept = mz_cal_params
        rdr = UReaderSetMzParams(uimf, _UREADER_WORKERS, mz_cal_slope, mz_cal_intercept)
    else:
        rdr = UReader(uimf, _UREADER_WORKERS)
    print('-' * 40)
    print('SURVEY ANALYSIS')
    print('-' * 40)
    # create a map between dataset identifiers and dataset labels
    dsid_to_dslbl = {
        row["dataset_id"]: row["dataset_label"] 
        for row in amt_cur.execute("SELECT dataset_id, dataset_label FROM Datasets").fetchall()
    }
    # load the targets from the AMT tag database
    target_labels, target_mzs = [], []
    qry_sel_targets = """--sqlite3
        SELECT 
            target_id, 
            target_labels, 
            target_mz, 
            peptide_ids, 
            dataset_ids
        FROM 
            Targets 
        WHERE 
            target_mz BETWEEN :min_tmz AND :max_tmz
        ORDER BY 
            target_mz
    ;"""
    for row in amt_cur.execute(qry_sel_targets, 
                               {"min_tmz": min_tmz, "max_tmz": max_tmz}
                               ).fetchall():
        # make sure the target is present in one of the selected datasets
        if _target_is_in_selected_dataset(sel_dsets, dsid_to_dslbl, row["dataset_ids"]):
            target_labels.append(f"{row["target_id"]}_{row["target_labels"]}")
            target_mzs.append(row["target_mz"])
    # analyze survey data to see which targets appear to be present
    survey_targets = _survey_analysis(rdr, 
                                      target_labels, 
                                      target_mzs, 
                                      mz_ppm, 
                                      frame_mod, 
                                      scan_mod, 
                                      i_min, 
                                      min_peak_ht, 
                                      min_peak_fwhm, 
                                      max_peak_fwhm, 
                                      min_dist)
    # filter the survey targets 
    survey_targets_keep = [[], [], [], []]
    for lbl, mz, dt, fwhm in zip(*survey_targets):
        add = True
        for klbl, kmz, kdt, kfwhm in zip(*survey_targets_keep):
            flag = lbl == klbl and abs(mz - kmz) <= tol_from_ppm(mz, mz_ppm) and abs(dt - kdt) <= min_dist
            if flag:
                add = False
                break
        if add:
            survey_targets_keep[0].append(lbl)
            survey_targets_keep[1].append(mz)
            survey_targets_keep[2].append(dt)
            survey_targets_keep[3].append(fwhm)
    # report how many survey targets there were
    print('*' * 40)
    print('survey targets (before filter):', len(survey_targets[0]))
    print('survey targets (after filter):', len(survey_targets_keep[0]))
    print('*' * 40)
    print('-' * 40)
    print('FULL ANALYSIS')
    print('-' * 40)
    # analyze full data near peaks found in survey data
    found_params = _full_analysis(rdr, 
                                  survey_targets_keep, 
                                  i_min, 
                                  mz_ppm, 
                                  min_peak_ht, 
                                  min_peak_fwhm, 
                                  max_peak_fwhm, 
                                  full_min_psnr, 
                                  True)
    print('-' * 40)
    return found_params


def targeted_extraction_no_amt(
        mbi_f: str,
        target_labels: List[str],
        target_mzs: List[float],
        mz_ppm: float, 
        scan_min: int,
        scan_max: int,
        min_peak_ht: float, 
        min_peak_fwhm: float, 
        max_peak_fwhm: float,
        min_psnr: float
    ):
    """
    """
    # initialize MOBIE reader
    rdr = MReader(mbi_f)
    # extract the sparse data
    sparse_data = rdr.extract_sparse_data(scan_min, scan_max)
    # search the data for the targets
    found = {}
    for lbl, mz in zip(target_labels, target_mzs):
        print("-" * 80)
        print(f"{lbl=} {mz=}")
        mz_tol = tol_from_ppm(mz, mz_ppm)
        atd = rdr.select_sparse_atd(sparse_data, mz, mz_tol, scan_min, scan_max)
        # look for up to 3 peaks
        pkdts, pkhts, pkwts = find_peaks_1d_gauss(*atd, 0.1, min_peak_ht, min_peak_fwhm, max_peak_fwhm, 3, True)
        peaks_found = False
        for pkdt, pkht, pkwt in zip(pkdts, pkhts, pkwts):
            psnr = calc_gauss_psnr(*atd, (pkdt, pkht, pkwt))
            if psnr >= min_psnr:
                peak_data = {
                    'mz': mz,
                    # convert scan to ms
                    'dt_ms': pkdt * rdr.avg_dt_bin_width,
                    'pk_params': (pkdt, pkht, pkwt, psnr),
                    'atd': atd,
                    'spectrum': rdr.select_sparse_spectrum(sparse_data, mz - 1.5, mz + 2.5, scan_min, scan_max)
                }
                if not peaks_found:
                    found[lbl] = [peak_data]
                else:
                    found[lbl].append(peak_data)
                peaks_found = True
                print(f"arrival time: {pkdt:.0f} +/- {pkwt:.0f} scans ({pkht:.3e}, {psnr:.1f})")
                # plot MS1 spectrum
                # ax = plot_spectrum(
                #     *peak_data['spectrum'],
                #     figsize=(4, 1.5),
                #     c="#0EB1B1"
                # )
                # ax.axvline(mz, c="#D82189", ls='--', lw=0.75, zorder=-1)
                # ax.axvline(mz + 1, c="#D82189", ls='--', lw=0.75, zorder=-1)
                # ax.axvline(mz + 2, c="#D82189", ls='--', lw=0.75, zorder=-1)
                # plt.show()
                # plt.close()
        if not peaks_found:
            print('no ATD peaks found')
        else:
            # plot atd
            # ax = plot_atd(
            #     *atd, 
            #     figsize=(3, 2), 
            #     dt_unit='scans',
            #     dt_range=(min(pkdts) - max(pkwts) * 4, max(pkdts) + max(pkwts) * 4),
            #     c="#0EB1B1"
            # )
            # for pkdt, pkht, pkwt in zip(pkdts, pkhts, pkwts):
            #     dt = atd[0]
            #     dtp = dt[(dt >= pkdt - 2 * pkwt) & (dt <= pkdt + 2 * pkwt)]
            #     ax.plot(dtp, _gauss(dtp, pkdt, pkht, pkwt), ls="-", c="#D82189", lw=1)
            # plt.show()
            # plt.close()
            pass
    print("-" * 80)
    # close the reader
    rdr.close()
    return found


def _check_isotopes(lbl: str, 
                    mz: float, 
                    peak_mzs: Tuple[float, float, float], 
                    peak_hts: Tuple[float, float, float]
                    ):
    """
    """
    unique_zs = set([int(_[-1]) for _ in lbl.split('|')])
    mzt = 0.1  # m/z tolerance for matching isotopologs
    evidence_for_zs = {}
    for z in unique_zs:
        # calculate putative isotopolog m/z for a given charge state
        mmz, m1mz, m2mz = mz, mz + 1 / z, mz + 2 / z
        # try to match them to observed peaks
        pmmz, pm1mz, pm2mz = None, None, None
        pmht, pm1ht, pm2ht = None, None, None
        for pmz, pht, in zip(peak_mzs, peak_hts):
            if abs(pmz - mmz) <= mzt:
                pmmz = pmz
                pmht = pht
            if abs(pmz - m1mz) <= mzt:
                pm1mz = pmz
                pm1ht = pht
            if abs(pmz - m2mz) <= mzt:
                pm2mz = pmz
                pm2ht = pht
        # decide if there is sufficient evidence for a charge state based on whether 
        # at least M and M+1 isotopologs were observed
        if pmmz is not None and pm1mz is not None:
            evidence_for_zs[z] = ((pmmz, pm1mz, pm2mz), (pmht, pm1ht, pm2ht))
    if len(evidence_for_zs) > 1:
        # if there is evidence for z=2 then there will be evidence for z=1 as well
        # since for z=2 M+2 is the same as M+1 for z=1
        # so if the two options are 1 and 2 then we can ignore z=1
        zs = [_ for _ in evidence_for_zs.keys()]
        if 3 not in zs:
            z = 2
            (pmmz, pm1mz, pm2mz), (pmht, pm1ht, pm2ht) = evidence_for_zs[z]
        else:
            return 0, "", (pmmz, None, None), (pmht, None, None)
            # ?
            #raise RuntimeError("should not have evidence for charge states (1 + 3) or (2 + 3)")
    elif len(evidence_for_zs) == 0:
        # insufficient evidence for any charge state in original labels
        return 0, "", (pmmz, None, None), (pmht, None, None) 
    # only evidence for 1 chanrge state, unpack it then filter the labels and return it
    else:
        z = [_ for _ in evidence_for_zs.keys()][0]
        (pmmz, pm1mz, pm2mz), (pmht, pm1ht, pm2ht) = evidence_for_zs[z]
    # check if this seems to be an isotopolog of something else, i.e. look for a peak
    # at M-1 for the selected charge state that is at least 90% the height of the 
    # putative M peak
    m_minus_1_mz, m_minus_1_ht = pmmz - 1 / z, None
    for pmz, pht in zip(peak_mzs, peak_hts):
        if abs(pmz - m_minus_1_mz) <= mzt:
            m_minus_1_ht = pht
    if m_minus_1_ht is not None and m_minus_1_ht > 0.9 * pmht:
        # this looks like an isotopolog of something else, return "" for the label
        return z, "",  (pmmz, pm1mz, pm2mz), (pmht, pm1ht, pm2ht)
    # filter labels based on evidence for charge states
    new_lbls = []
    for seq_z in lbl.split("|"):
        seq, _z, = seq_z[:-2], seq_z[-1:]
        seq = seq.split("_")[1] if "_" in seq else seq
        if int(_z) == z:
            new_lbls.append(seq)
    new_lbl = "|".join(new_lbls)
    return z, new_lbl, (pmmz, pm1mz, pm2mz), (pmht, pm1ht, pm2ht)


def post_process_results(results, out_file):
    """ 
    post-process feature extraction results by checking isotope 
    distributions and write filtered IDs to CSV 
    """
    print('-' * 40)
    print('POST-PROCESS RESULTS')
    print('-' * 40)
    s = "{},{},{:.4f},{:.1f},{:.1f},{:.4e},{:.1f}\n"
    with open(out_file, "w") as out:
        out.write("label,z,mz,arrival_time,fwhm,height,psnr\n")
        for lbl in results.keys():
            print("- " * 20)
            print(lbl)
            for i in range(len(results[lbl])):
                #print("peak", i + 1)
                #atd = results[lbl][i]["atd"]
                # # plot atd
                # ax = plot_atd(*atd, figsize=(6, 2))
                dt, ht, wt, snr = results[lbl][i]["pk_params"]
                # ax.axvline(dt, ls="--", c='k', lw=0.75)
                # ax.set_xlim([dt - 1000, dt + 1000])
                # plt.show()
                # plt.close()
                spectrum = results[lbl][i]["spectrum"]
                mz = results[lbl][i]["mz"]
                # find peaks in mass spectrum
                pkmzs, pkhts, _ = find_peaks_1d_localmax(*spectrum, 0.1, 1e6, 0.05, 0.3, 0.25)
                # check isotopes
                iso_z, iso_lbl, *_ = _check_isotopes(lbl, results[lbl][i]["mz"], pkmzs, pkhts)
                print("new label:", iso_lbl)
                if iso_lbl != "":
                    # only write results that pass through the check_isotopes function
                    # indicated by the returned label not being an empty string
                    out.write(s.format(iso_lbl, iso_z, mz, dt, wt, ht, snr))
                # # plot spectrum
                # ax = plot_spectrum(*spectrum, figsize=(6, 2))
                # ax.plot(pkmzs, pkhts, 'kx', ms=7, mew=1.5)
                # for imz, iht in zip(iso_pkmzs, iso_pkhts):
                #     if imz is not None:
                #         ax.plot(imz, iht, 'go', ms=7, mew=1.5, fillstyle='none')
                # ax.axvline(mz, ls="--", c='k', lw=0.75)
                # plt.show()
                # plt.close()
