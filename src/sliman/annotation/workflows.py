"""
sliman/annotation/workflows.py

Dylan Ross (dylan.ross@pnnl.gov)

    High-level workflow functions for feature annotation
"""


from itertools import product, combinations

import numpy as np

from mzapy.isotopes import ms_adduct_formula, predict_m_m1_m2

from sliman.annotation.peptides import Peptide


def _consolidate_matches_2(m_feats, m1_feats, am, am1, am2, dt_tol, dot_threshold):
    """ 
    go through all combinations of potential matches (M, M+1, M+2), yield any valid combinations
    
    valid combinations are ones that satisfy the following conditions:
    * arrival times must all be within dt_tol of one another
    * relative abundances of isotopologues must be similar to predicted (dot product <= dot_threshold)
    """
    for m_feat, m1_feat in product(m_feats, m1_feats):
        # maximum difference in DT between this group of features
        max_ddt = max([abs(dt_a - dt_b) for dt_a, dt_b in combinations([m_feat[1], m1_feat[1]], 2)])
        print('max dDT: {:.0f}'.format(max_ddt))
        if max_ddt <= dt_tol:
            # check relative abundances
            f_abs = [m_feat[2], m1_feat[2], 0]
            f_ras = np.array(f_abs) / np.linalg.norm(f_abs)
            p_abs = [am, am1, am2]
            p_ras = np.array(p_abs) / np.linalg.norm(p_abs)
            dp = np.dot(f_ras, p_ras)
            print('obs  RA: [{:.3f}, {:.3f}]'.format(*f_ras))
            print('pred RA: [{:.3f}, {:.3f}]'.format(*p_ras))
            print('dot product: {:.3f}'.format(dp))
            if dp >= dot_threshold:
                yield {'M': m_feat, 'M+1': m1_feat}


def _consolidate_matches_3(m_feats, m1_feats, m2_feats, am, am1, am2, dt_tol, dot_threshold):
    """ 
    go through all combinations of potential matches (M, M+1, M+2), yield any valid combinations
    
    valid combinations are ones that satisfy the following conditions:
    * arrival times must all be within dt_tol of one another
    * relative abundances of isotopologues must be similar to predicted (dot product <= dot_threshold)
    """
    for m_feat, m1_feat, m2_feat in product(m_feats, m1_feats, m2_feats):
        # maximum difference in DT between this group of features
        max_ddt = max([abs(dt_a - dt_b) for dt_a, dt_b in combinations([m_feat[1], m1_feat[1], m2_feat[1]], 2)])
        print('max dDT: {:.0f}'.format(max_ddt))
        if max_ddt <= dt_tol:
            # check relative abundances
            f_abs = [m_feat[2], m1_feat[2], m2_feat[2]]
            f_ras = np.array(f_abs) / np.linalg.norm(f_abs)
            p_abs = [am, am1, am2]
            p_ras = np.array(p_abs) / np.linalg.norm(p_abs)
            dp = np.dot(f_ras, p_ras)
            print('obs  RA: [{:.3f}, {:.3f}, {:.3f}]'.format(*f_ras))
            print('pred RA: [{:.3f}, {:.3f}, {:.3f}]'.format(*p_ras))
            print('dot product: {:.3f}'.format(dp))
            if dp >= dot_threshold:
                yield {'M': m_feat, 'M+1': m1_feat, 'M+2': m2_feat}
            

def match_peptide_iso_dist(seq, z, feats, mz_tol, dt_tol, dot_threshold, require_m2=False, m2_lower_limit=0.1):
    """ 
    find a set of features that match the expected isotope distribution of a peptide 
    
    Parameters
    ----------
    seq : ``str``
        peptide sequence
    z : ``int``
        peptide charge state
    feats : ``list(...)``
        unlabeled 2D features, each row has m/z, DT, abundance, and the rest of the attributes are ignored
    mz_tol : ``float``
        m/z tolerance for matching features
    dt_tol : ``float``
        DT tolerance for matching features
    dot_threshold : ``float``
        minimum dot product for matching isotope distribution
    require_m2 : ``bool``, defaut=False
        require M, M+1, M+2 isotopologues to be found (if False, M+2 can be ignored if abundance is below m2_lower_limit)
    m2_lower_limit : ``float``, default=0.1
        minimum absolute abundance for M+2 isotopologue to be able to ignore it

    Yields
    ------
    consolidated_feature : ``dict(...)``
        dictionary with matched features for M, M+1, and M+2 isotopologues
    """
    if z > 0:
        adduct = '[M+{z}H]{z}+'.format(z=z) if z > 1 else '[M+H]+'
    else:
        adduct = '[M-{z}H]{z}-'.format(z=-z) if z < -1 else '[M-H]-'
    pep = Peptide(seq)
    formula = ms_adduct_formula(pep.formula, adduct)
    masses, abuns = predict_m_m1_m2(formula, relative_abundance=False)
    mz, mz1, mz2 = [m / z for m in masses]
    am, am1, am2 = abuns
    # find each isotopologue
    # M
    m_feats = []
    for fmz, fdt, fab, *_ in feats:
        if abs(mz - fmz) <= mz_tol:
            print('{:>20s} z={:1d} {:3s} {:8.3f} {:.3f} {:8.3f} {:5.0f} {:5.3e}'.format(seq, z, 'M', mz, am, fmz, fdt, fab))
            m_feats.append([fmz, fdt, fab])
    if m_feats != []:  # M isotopologue found 
        # M+1
        m1_feats = []
        for fmz, fdt, fab, *_ in feats:
            if abs(mz1 - fmz) <= mz_tol:
                print('{:>20s} z={:1d} {:3s} {:8.3f} {:.3f} {:8.3f} {:5.0f} {:5.3e}'.format(seq, z, 'M+1', mz1, am1, fmz, fdt, fab))
                m1_feats.append([fmz, fdt, fab])
        if m1_feats != []:  # M+1 isotopologue found
            m2_feats = []
            for fmz, fdt, fab, *_ in feats:
                if abs(mz2 - fmz) <= mz_tol:
                    print('{:>20s} z={:1d} {:3s} {:8.3f} {:.3f} {:8.3f} {:5.0f} {:5.3e}'.format(seq, z, 'M+2', mz2, am2, fmz, fdt, fab))
                    m2_feats.append([fmz, fdt, fab])
            if m2_feats != []:
                for consolidated in _consolidate_matches_3(m_feats, m1_feats, m2_feats, am, am1, am2, dt_tol, dot_threshold):
                    yield consolidated
            elif (not require_m2) and am2 <= m2_lower_limit:
                for consolidated in _consolidate_matches_2(m_feats, m1_feats, am, am1, am2, dt_tol, dot_threshold):
                    yield consolidated

