"""
sliman/_util.py

Dylan Ross (dylan.ross@pnnl.gov)

    Internal module with general utilities
"""


def tol_from_ppm(mz: float, ppm: float) -> float :
    return mz * ppm / 1e6

