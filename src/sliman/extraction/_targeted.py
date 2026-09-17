"""
sliman/extraction/targeted.py

Dylan Ross (dylan.ross@pnnl.gov)

    Internal module with utilities for targeted analyses
"""


import numpy as np


def select_survey_atd_sparse(target_mz, mz_tol, n_scans, scan_mod, survey_data_sparse):
    """ 
    select out ATD for m/z from sparse survey data, sum together using specified scan_mod
    
    Parameters
    ----------
    target_mz : ``float``
        target m/z
    mz_tol : ``float``
        m/z tolerance
    n_scans : ``int``
        number of scans in a frame
    scan_mod : ``int``
        determines how densely the arrival time (scans) is sampled in survey data, 
        must match the value used to extract survey data
    survey_data_sparse : ``numpy.ndarray(...)``
        survey data, 3 columns: scan, m/z, intensity
    
    Returns
    -------
    atd_scans : ``numpy.ndarray(int)``
    atd_intensities : ``numpy.ndarray(float)``
        m/z-selected ATD
    """
    ss = np.arange(0, n_scans, scan_mod)
    si = np.zeros(ss.shape)
    idx = (survey_data_sparse[1] >= target_mz - mz_tol) & (survey_data_sparse[1] <= target_mz + mz_tol)
    for s, i in zip(survey_data_sparse[0][idx], survey_data_sparse[2][idx]):
        si[int(s) // scan_mod] += i
    return ss, si 


def select_full_atd_sparse(target_mz, mz_tol, n_scans, full_data_sparse):
    """ 
    select out ATD for m/z from sparse full data, sum together
    
    Parameters
    ----------
    target_mz : ``float``
        target m/z
    mz_tol : ``float``
        m/z tolerance
    n_scans : ``int``
        number of scans in a frame
    full_data_sparse : ``numpy.ndarray(...)``
        full data, 3 columns: scan, m/z, intensity
    
    Returns
    -------
    atd_scans : ``numpy.ndarray(int)``
    atd_intensities : ``numpy.ndarray(float)``
        m/z-selected ATD
    """
    # the function logic is the same as _select_survey_atd_sparse with 
    # scan_mod set to 1 (i.e., full scan range), just use that
    return select_survey_atd_sparse(target_mz, mz_tol, n_scans, 1, full_data_sparse)


def select_full_spectrum_sparse(target_mz, mz_bin_size, scan_min, scan_max, full_data_sparse):
    """
    select out mass spectrum for arrival time range around target m/z from sparse full data
    m/z range is target - 1.5 to target + 2.5, this covers main isotopes for all charge states

    Parameters
    ----------
    target_mz : ``float``
        target m/z
    mz_bin_size : ``float``
        size of m/z for summing together spectrum
    scan_min : ``int``
    scan_max : ``int``
        scan bounds for arrival time selection
    full_data_sparse : ``numpy.ndarray(...)``
        full data, 3 columns: scan, m/z, intensity

    Returns
    -------
    ms1_mz : ``numpy.ndarray(float)``
    ms1_intensities : ``numpy.ndarray(float)``
        arrival time-selected mass spectrum
    """
    mz_min, mz_max = target_mz - 1.5, target_mz + 2.5
    sm = np.arange(mz_min, mz_max + mz_bin_size, mz_bin_size)
    si = np.zeros(sm.shape)
    m_to_idx = lambda m: int((m - mz_min) // mz_bin_size)
    idx = (full_data_sparse[1] >= mz_min) & (full_data_sparse[1] <= mz_max) & \
          (full_data_sparse[0] >= scan_min) & (full_data_sparse[0] <= scan_max)
    for m, i in zip(full_data_sparse[1][idx], full_data_sparse[2][idx]):
        si[m_to_idx(m)] += i
    return sm, si 

