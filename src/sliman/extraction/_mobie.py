"""
sliman/extraction/_mobie.py
Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for extracting from MOBIE data
"""


from typing import Any, List
from time import time

import numpy as np
from scipy.sparse import coo_array

from mbisdk import *


class MReader:
    """
    Attributes
    ----------
    mbi_f : `str`
        .mbi data file
    mbi : `MBIFile`
        Mobie SDK interface object 
    """

    def __init__(self, mbi_f: str):
        """ 
        Parameters
        ----------
        mbi_f
            .mbi data file
        """
        # store attributes
        self.__mbi_f = mbi_f
        # connect to the MBIFile
        self.__mbi = MBIFile(self.mbi_f)
        self.mbi.Init()
        if not self.mbi.IsInitialized():
            raise RuntimeError(f"unable to initialize MBIFile from {self.mbi_f}")
        # set flag to indicate the instance is ready to read data
        self.__ready = True
        # store some useful data stats
        self.__n_frames = self.mbi.GetNumFrames()
        self.__scans_per_frame = list(map(self.mbi.GetMaxScansInFrame,  range(1, self.n_frames + 1)))
        self.__max_points_in_scan = self.mbi.GetMaxPointsInScan()
        self.__all_mzs = np.array(list(map(self.mbi.GetCalibration().IndexToMz, range(self.max_points_in_scan))))
        self.__avg_at_bin_width = np.mean([
            self.mbi.GetFrameMetadata(f_num).ReadDouble("frm-dt-period") 
            for f_num in range(1, self.n_frames + 1)
        ])

    # --- attribute getters ---

    @property
    def mbi_f(self) -> str : return self.__mbi_f

    @property
    def mbi(self) : return self.__mbi

    @property
    def ready(self) -> bool : return self.__ready

    @property
    def n_frames(self) -> int : return self.__n_frames

    @property
    def scans_per_frame(self) -> List[int] : return self.__scans_per_frame

    @property
    def max_points_in_scan(self) -> int : return self.__max_points_in_scan

    @property
    def all_mzs(self) -> Any : return self.__all_mzs

    @property
    def avg_dt_bin_width(self) -> float : return self.__avg_at_bin_width
        
    # --- normal methods ---

    def extract_sparse_data(self, scan_min: int, scan_max: int):
        """
        """
        if not self.ready:
            raise RuntimeError("this instance is not ready to read data from the MOBIE file")
        print("reading scan data ... ")
        t0 = time()
        scans, mzbins, iis = [], [], []
        n = 0
        for i_frame in range(self.n_frames):
            print(f"\r\tframe: {i_frame + 1}", end="")
            frame = self.mbi.GetFrame(i_frame + 1)
            _scan_max = min(scan_max + 1, self.scans_per_frame[i_frame])
            for scan in range(scan_min, _scan_max + 1):
                for mzidx, intensity in (
                    lambda spec: zip(spec.indices, spec.intensities)
                )(frame.GetMassSpectrum(scan)):
                    scans.append(scan)
                    mzbins.append(mzidx)
                    iis.append(intensity)
                    n += 1
            frame.Unload()
        print(f"\n{n=}")
        print("converting to CSR array ... ")
        sparse_data = coo_array((iis, (mzbins, scans))).tocsr()
        print(f"... done (elapsed: {time() - t0:.3f} s)")
        return sparse_data
    
    def select_sparse_atd(self, sparse_data, mz: float, tol: float, scan_min: int, scan_max: int):
        """
        """
        mz_idx_min = np.argmin(np.abs(self.all_mzs - (mz - tol)))
        mz_idx_max = np.argmin(np.abs(self.all_mzs - (mz + tol)))
        intensities = np.sum(sparse_data[mz_idx_min:mz_idx_max, :].todense(), axis=0)[scan_min:scan_max + 1]
        scans = np.arange(scan_min, scan_max + 1)
        # ensure that ATD always has the same number of points in scans and intensities arrays
        if (padding_needed := len(scans) - len(intensities)) > 0:
            intensities = np.pad(intensities, (0, padding_needed), "constant")
        return scans, intensities
    
    def select_sparse_spectrum(self, sparse_data, mz_min: float, mz_max: float, scan_min: int, scan_max: int):
        """
        """
        mz_idx_min = np.argmin(np.abs(self.all_mzs - mz_min))
        mz_idx_max = np.argmin(np.abs(self.all_mzs - mz_max))
        intensities = np.sum(sparse_data[mz_idx_min:mz_idx_max, scan_min:scan_max].todense(), axis=1)
        return self.all_mzs[mz_idx_min:mz_idx_max], intensities
    
    def close(self):
        """ close the underlying `MBIFile` interface, unset `self.ready` flag """
        self.__ready = False
        self.mbi.Close()
        