"""
sliman/annotation/calibration/__init__.py

Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for performing reference-free CCS calibration and filtering of peptide annotations

    NOTE: The code is here in the __init__.py module because it was originally all contained within
          a single module, but it ended up getting too big so I turned it into a package and moved
          some stuff into internal modules within that package. This setup makes it so functions
          can be imported like ``from sliman.annotation.calibration import foo, bar`` insead of
          having to refer to a submodule like 
          ``from sliman.annotation.calibration.something import foo, bar``
"""


import sqlite3
import pickle
import os
from typing import Callable, Dict, List

import numpy as np
from numpy import typing as npt
from mzapy.calibration import CCSCalibrationTW
from mzapy.view import plot_tw_ccs_calibration

from sliman.annotation.calibration import _queries
from sliman.annotation.calibration import _plots


def _insert_log_entry(dbf: str, entry: str) -> None :
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    cur.execute(_queries.QRY_INS_LOG_ENTRY, (entry,))
    con.commit()
    con.close()


def _init_peptide_db(dbf: str, overwrite: bool) -> None :
    """
    initialize the peptide database, possibly overwriting it
    """
    # check if the DB file exists, possibly overwrite
    if os.path.isfile(dbf):
        if overwrite:
            os.remove(dbf)
        else:
            msg = f"_init_peptide_db: DB file {dbf} already exists (and overwrite=False)"
            raise RuntimeError(msg)
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    cur.executescript(_queries.PEPTIDE_DB_SCHEMA)
    # cleanup
    con.commit()
    con.close()
    # add an entry to the log
    _insert_log_entry(dbf, "_init_peptide_db")


def _add_results_to_peptide_db(dbf: str, results_file: str, sample_name: str, 
                               scan_to_ms_cb: Callable[[float], float]
                               ) -> None :
    """
    add a single set of results into the peptide database 
    """
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # add peptides from results file into DB
    # keep track of peptides/charge states that have already been added
    # map (sequence, z) to rowid in Peptides table
    seq_z_ids = {(seq, z): pep_id for pep_id, seq, z in cur.execute(_queries.QRY_SEL_PEPS_1)}
    with open(results_file, "r") as f:
        _ = next(f)
        for line in f.readlines():
            # parse the line, make necessary type conversions
            lbl, z, mz, arrival_time_scan, fwhm, height, psnr = line.split(",")
            z = int(z)
            mz, arrival_time_scan, fwhm, height, psnr = [float(_) for _ in [mz, arrival_time_scan, fwhm, height, psnr]]
            # add an entry to the Features table
            qdata = (None, sample_name, mz, scan_to_ms_cb(arrival_time_scan), fwhm, height, psnr)
            cur.execute(_queries.QRY_INSERT_FEAT, qdata)
            feat_id = cur.lastrowid
            # add entries to the Peptides table
            for seq in lbl.split("|"):
                pep_id = seq_z_ids.get((seq, z))
                if pep_id is None:
                    # add a new entry
                    cur.execute(_queries.QRY_INSERT_PEP, (None, seq, z))
                    pep_id = cur.lastrowid
                    seq_z_ids[(seq, z)] = pep_id
                # pep_id now has a valid ID, whether from the cache or 
                # a newly added entry
                cur.execute(_queries.QRY_INSERT_FEATPEPMAP, (feat_id, pep_id))
    # cleanup
    con.commit()
    con.close()
    # add an entry to the log
    _insert_log_entry(dbf, f"_add_results_to_peptide_db: {sample_name=}")


def setup_peptide_db(dbf: str, results: Dict[str, str], 
                     scan_to_ms_cb: Callable[[float], float], 
                     overwrite: bool = False
                     ) -> None :
    """
    Initialize peptide database and fill with putative peptide IDs from 
    post-processed analysis results

    Parameters
    ----------
    dbf : ``str``
        path to database file
    results : ``dict(str:str)``
        dict of post-processed results to include, mapping sample name to corresponding 
        results csv file name
    scan_to_ms_cb : ``func(float) -> float``
        function for mapping arrival time scans into milliseconds
    overwrite : ``bool``, default=False
        if the specified database file already exists prior to calling this function, it
        will be overwritten if this flag is set to True, otherwise a RuntimeError will
        be raised
    """
    # init the database
    _init_peptide_db(dbf, overwrite)
    # iterate through datasets and add them to database
    n = 0
    for sample_name, results_file in results.items():
        _add_results_to_peptide_db(dbf, results_file, sample_name, scan_to_ms_cb)
        n += 1
    # add an entry to the log
    _insert_log_entry(dbf, f"setup_peptide_db: {n=} analysis results files")


def gen_dl_pred_input(dbf: str, dl_input_file: str, omit_symbols: bool = True) -> None :
    """
    fetch peptide data from the database and construct input file for DL 
    CCS prediction

    .. note::
        When the ``omit_symbols`` flag is set, the symbols *, #, and @ which commonly 
        follow residues M, R, and K, respectively, as the DL model cannot handle these. 
        This will impart a slight mismatch between the sequence represented in the DL 
        input (and output) and what is stored in the peptide database, but only the 
        peptide ID is used to assign the predicted values so that is ok. The other thing 
        is that these symbols do represent chemical modifications to the peptide at the 
        specified residues, therefore omitting them will probably introduce some amount 
        of error in the predicted CCS value from the DL model. We are assuming that this 
        effect is very small, and the value of having the predicted CCS values for these 
        entries outweighs the drawbacks of potentially introducing small errors in the 
        predicted CCS values from omission of these modifications.

    Parameters
    ----------
    dbf : ``str``
        path to database file
    dl_input_file : ``str``
        path to save the DL input file (a CSV)
    omit_symbols : ``bool``, default=True
        if True, omit the symbols *, #, and @ which commonly follow residues M, R, and K, 
        respectively, as the DL model cannot handle these.
    """
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # fetch the peptides and charge states from the database
    # and write to .csv
    with open(dl_input_file, "w") as f:
        f.write("pep_id,Modified sequence,Charge\n")
        for pep_id, seq, z in cur.execute(_queries.QRY_SEL_PEPS_1):
            _seq = seq.split(".")[1]
            if omit_symbols:
                _seq = _seq.replace("*", "").replace("#", "").replace("@", "")
            f.write(f"{pep_id},_{_seq}_,{z}\n")
    # cleanup
    con.commit()
    con.close()
    # add an entry to the log
    _insert_log_entry(dbf, f"gen_dl_pred_input: {dl_input_file=}")


def add_dl_pred_ccs_to_db(dbf: str, dl_output_file: str) -> None :
    """
    Add DL predicted CCS values to database
    """
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # load the candidates from the database for each charge state
    # and write to .csv
    with open(dl_output_file, "r") as f:
        _ = next(f)
        for line in f.readlines():
            # Should we rely on the specific structure of the CSV here? Yes.
            _, pep_id, *_, ccs = line.strip().split(",")
            cur.execute(_queries.QRY_INS_PEP_DL_CCS, (int(pep_id), float(ccs)))
    # cleanup
    con.commit()
    con.close()
    # add an entry to the log
    _insert_log_entry(dbf, f"add_dl_pred_ccs_to_db: {dl_output_file=}")


def _sample_calibrants(dts: npt.ArrayLike) -> npt.ArrayLike:
    """ 
    Return indices of randomly sampled calibrants from an array of arrival times.
    Sampling is performed such that calibrants from across the arrival time range 
    are represented. 
    """
    sampled = []
    n = len(dts)
    indices = np.arange(n)
    quintiles_ish = np.percentile(dts, [10, 25, 40, 60, 75, 90])
    for i in range(5):
        idx = (dts >= quintiles_ish[i]) & (dts <= quintiles_ish[i + 1])
        # select a min of 2, max of 0.2 * n calibrants per quintile
        n = max(int(len(indices[idx]) * 0.2), 2)
        sampled += [_ for _ in np.random.choice(indices[idx], n, replace=False)]
    return sampled


def _store_calibrant(cur: sqlite3.Cursor, 
                     z: int, 
                     seq: str, 
                     mz: float, 
                     dt: float, 
                     ccs: float
                     ) -> None :
    """
    store calibrants for a single charge state in the database
    """
    cur.execute(_queries.QRY_INS_CALIBRANT, 
                (z, seq, mz, dt, ccs))


def _select_calibrants(dbf: str, zs: List[int], seed: int) -> None :
    np.random.seed(seed)
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # load the candidates from the database for each charge state
    seqs = {}
    mzs = {}
    dts = {}
    ccss = {}
    mzs_all = {}
    dts_all = {}
    ccss_all = {}
    for z in zs:
        for d in [seqs, mzs, dts, ccss, mzs_all, dts_all, ccss_all]:
            d[z] = []
        for seq, mz, dt, ccs in cur.execute(_queries.QRY_SEL_CANDIDATES, (z,)):
            seqs[z].append(seq)
            mzs[z].append(mz)
            dts[z].append(dt)
            ccss[z].append(ccs)
        for seq, mz, dt, ccs in cur.execute(_queries.QRY_SEL_ALL_FEATS, (z,)):
            mzs_all[z].append(mz)
            dts_all[z].append(dt)
            ccss_all[z].append(ccs)
    # plot the potential calibrant points (mz, dt, ccs)
    #_plots.plot_candidate_mz_at_ccs(mzs, dts, ccss)
    #_plots.plot_candidate_at_hists(dts)
    # select the calibrants
    sampled = {}
    for z in zs:
        #print("=" * 15, f"z={z}", "=" * 15)
        #print("calibrants")
        sampled[z] = _sample_calibrants(dts[z])
        for idx in sampled[z]:
            seq, mz, dt, ccs = seqs[z][idx], mzs[z][idx], dts[z][idx], ccss[z][idx]
            #print(f"{seqs[z][idx]:>40s} {mzs[z][idx]:9.4f} {dts[z][idx]:.3f} {ccss[z][idx]:.2f}")
            _store_calibrant(cur, z, seq, mz, dt, ccs)
    # plot sampled calibrants
    #_plots.plot_calibrant_mz_at_ccs(mzs_all, dts_all, ccss_all, mzs, dts, ccss, sampled)
    # cleanup
    con.commit()
    con.close()
    # add a log entry
    _insert_log_entry(dbf, "_select_calibrants")


def _create_calibrations(dbf: str, zs: List[int]) -> None :
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    mzs, dts, ccss = {}, {}, {}
    cals = {}
    for z in zs:
        _mzs = []
        _dts = []
        _ccss = []
        for mz, dt, ccs in cur.execute(_queries.QRY_SEL_CALIBRANTS, (z,)):
                _mzs.append(float(mz))
                _dts.append(float(dt))
                _ccss.append(float(ccs))
        mzsa, dtsa, ccssa = np.array([_mzs, _dts, _ccss]) 
        idx = np.argsort(dtsa)
        cal = CCSCalibrationTW(mzsa[idx], dtsa[idx], ccssa[idx], z, "linear")
        #plot_tw_ccs_calibration(cal, figname=f"_figures/calibration_curve_z{z}.png")
        mzs[z], dts[z], ccss[z] = mzsa, dtsa, ccssa
        cals[z] = cal
        cur.execute(_queries.QRY_INS_CALIBRATION, (z, pickle.dumps(cal)))
    #_plots.plot_calibrant_at_ccs(mzs, dts, ccss, cals)
    # cleanup
    con.commit()
    con.close()
    # add a log entry
    _insert_log_entry(dbf, "_create_calibrations")


def _apply_calibrations(dbf: str, zs: List[int]):
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # load calibrations
    cals = {}
    for z, pkl in cur.execute(_queries.QRY_SEL_CALIBRATIONS):
        cals[z] = pickle.loads(pkl)
    # load the candidates from the database for each charge state
    for feat_id, mz, dt, z in cur.execute(_queries.QRY_SEL_FEATS).fetchall():
        ccs = cals[z].calibrated_ccs(mz, dt)
        cur.execute(_queries.QRY_INS_CCS, (feat_id, ccs))
    # cleanup
    con.commit()
    con.close()
    # add a log entry
    _insert_log_entry(dbf, "_apply_calibrations")


def _clear_previous_calibrations(dbf: str):
    """
    clear out any existing previous calibrations by deleting entries from
    Calibrants, Calibrations, and CalibratedCCS tables
    """
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # clear out the tables associated with CCS calibration
    cur.execute(_queries.QRY_RESET_CALIBRANTS)
    cur.execute(_queries.QRY_RESET_CALIBRATIONS)
    cur.execute(_queries.QRY_RESET_CALIBRATED_CCS)
    # cleanup
    con.commit()
    con.close()
    # add a log entry
    _insert_log_entry(dbf, "_clear_previous_calibrations")


def create_and_apply_ccs_calibrations(dbf: str, zs: List[int], seed: int = 42069) -> None :
    """
    Create CCS calibrations for each charge state present in the peptide
    database, then use them to generate calibrated CCS values. Intermediate
    information like the selected calibrants and individual calibration objects 
    are all stored within the database along with the calibrated CCS values. 

    Parameters
    ----------
    dbf : ``str``
        path to database file
    """
    # overwrite any existing calibrations if present
    _clear_previous_calibrations(dbf)
    # select calibrants from features in the database
    _select_calibrants(dbf, zs, seed)
    # create the calibrations (one for each charge state)
    _create_calibrations(dbf, zs)
    # apply the calibrations to calculate CCS for features in the database
    _apply_calibrations(dbf, zs)
    # add a log entry
    _insert_log_entry(dbf, f"create_and_apply_ccs_calibrations: {zs=} {seed=}")


def write_filtered_results(dbf, filtered_results_file, ccs_percent_threshold):
    """
    write rfCCS filtered peptide IDs to CSV

    Parameters
    ----------
    dbf : ``str``
        path to database file
    """
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # write everything to file
    hdr = "sample,peps,mz,z,at,cal_ccs,pred_ccs,pkht,n\n"
    line_template = "{smpl},{peps},{mz:.4f},{z},{at:.3f},{cal_ccs:.2f},{pccs},{pkht:.3e},{n}\n"
    filtered_results_file_n1 = filtered_results_file[:-4] + "_n1.csv"
    filtered_results_file_n2 = filtered_results_file[:-4] + "_n2.csv"
    with open(filtered_results_file, "w") as out, open(filtered_results_file_n1, "w") as out_n1, \
            open(filtered_results_file_n2, "w") as out_n2:
        out.write(hdr)
        out_n1.write(hdr)
        out_n2.write(hdr)
        # iterate through features first
        for feat_id, smpl, mz, dt, pkht, cal_ccs in cur.execute(_queries.QRY_SEL_FEATS_2).fetchall():
            # next iterate through peptides
            peps, zs, pccss = [], [], []
            for seq, z, pred_ccs in cur.execute(_queries.QRY_SEL_PEPS_2, (feat_id,)).fetchall():
                if abs(cal_ccs - pred_ccs) / pred_ccs <= ccs_percent_threshold / 100.:
                    peps.append(seq)
                    zs.append(z)
                    pccss.append(pred_ccs)
            assert len(set(zs)) <= 1, "post-processed features should only map to a single charge state (or none)"
            data = {
                "smpl": smpl, 
                "peps": "|".join(peps),
                "mz": mz,
                "z": zs[0] if len(zs) > 0 else None,
                "at": dt,
                "cal_ccs": cal_ccs,
                "pccs": "|".join([f"{_:.2f}" for _ in pccss]),
                "pkht": pkht,
                "n": len(peps)
            }
            if data["n"] > 0:
                out.write(line_template.format(**data))
                if data["n"] == 1:
                    out_n1.write(line_template.format(**data))
                if data["n"] <= 2:
                    out_n2.write(line_template.format(**data))
    # cleanup
    con.close()


def write_unfiltered_results(dbf, unfiltered_results_file):
    """
    """
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # write everything to file
    hdr = "sample,peps,mz,z,at,cal_ccs,pred_ccs,pkht,n\n"
    line_template = "{smpl},{peps},{mz:.4f},{z},{at:.3f},{cal_ccs:.2f},{pccs},{pkht:.3e},{n}\n"
    with open(unfiltered_results_file, "w") as out:
        out.write(hdr)
        # iterate through features first
        for feat_id, smpl, mz, dt, pkht, cal_ccs in cur.execute(_queries.QRY_SEL_FEATS_2).fetchall():
            # next iterate through peptides
            peps, zs, pccss = [], [], []
            for seq, z, pred_ccs in cur.execute(_queries.QRY_SEL_PEPS_2, (feat_id,)).fetchall():
                peps.append(seq)
                zs.append(z)
                pccss.append(pred_ccs)
            assert len(set(zs)) == 1, "post-processed features should only map to a single charge state"
            data = {
                "smpl": smpl, 
                "peps": "|".join(peps),
                "mz": mz,
                "z": zs[0],
                "at": dt,
                "cal_ccs": cal_ccs,
                "pccs": "|".join([f"{_:.2f}" for _ in pccss]),
                "pkht": pkht,
                "n": len(peps)
            }
            out.write(line_template.format(**data))
    # cleanup
    con.close()


def plot_filtered_peps(dbf, ccs_percent_threshold, figname):
    """
    """
    # connect to db
    con = sqlite3.connect(dbf)
    cur = con.cursor()
    # accumulate data 
    combined = {}
    # store feature arrival times and ccss
    dts_all = {}
    ccss_all = {}
    dts_filtered = {}
    ccss_filtered = {}
    info_filtered = {}
    # iterate through features first
    for feat_id, smpl, mz, dt, pkht, cal_ccs in cur.execute(_queries.QRY_SEL_FEATS_2).fetchall():
        # get the peptide IDs
        peps = []
        pccs = []
        for seq, z, pred_ccs in cur.execute(_queries.QRY_SEL_PEPS_2, (feat_id,)).fetchall():
            for d in [dts_all, ccss_all, dts_filtered, ccss_filtered, info_filtered]:
                if d.get(z) is None:
                    d[z] = []
            dts_all[z].append(dt)
            ccss_all[z].append(pred_ccs)
            if abs(cal_ccs - pred_ccs) / pred_ccs <= ccs_percent_threshold / 100.:
                peps.append(seq)
                pccs.append(pred_ccs)
        n = len(peps)
        peps = "|".join(peps)
        peps_pccs = "|".join([f"{_:.2f}" for _ in pccs])
        if len(peps) > 0:
            k = (mz, z)
            combined[k] = [smpl, peps, peps_pccs, dt, cal_ccs, pkht, n]
            if n == 1:
                # this introduces duplicates, acutally doesnt change the plots at all
                # because it ends up replotting over some of the same points so no need
                # to change but just FYI
                info_filtered[z].append(seq)
                dts_filtered[z].append(dt)
                ccss_filtered[z].append(pccs[0])
    # plot the filtered features
    #_plots.plot_filtered_at_ccs(dts_all, ccss_all, dts_filtered, ccss_filtered, figname, z_colors=False)
    _plots.plot_filtered_at_ccs(dts_all, ccss_all, dts_filtered, ccss_filtered, figname)

