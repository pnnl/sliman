"""
sliman/extraction/_uimf.py

Dylan Ross (dylan.ross@pnnl.gov)

    Internal module with minimal utilities for reading data from UIMF files
"""


from sqlite3 import connect
from multiprocessing import JoinableQueue, Process
import queue
import os
from typing import Optional, Set

import lzf
import numpy as np


# NOTE (Dylan Ross): This prevents deadlocks that occur when too much data 
#                    gets put into a queue before items are removed. It is 
#                    NOT a problem of the queue's max_size parameter (i.e. 
#                    the specified capacity of the queue in terms of item 
#                    count), which theoretically should be infinite when 
#                    set to 0. Instead, this problem arrises from the limited 
#                    capacity of the underlying pipe used to feed data between
#                    processes. For items with larger sizes (bytes), the pipe 
#                    can get filled up by relatively low item counts. This ulimately 
#                    causes a complete deadlock that prevents new items from being 
#                    added to the queue. Feeding chunks of spectra to be 
#                    unpacked should circumvent this issues by ensuring that 
#                    there are never too many items in either queue at one time
# batch size for unpacking (decompressing and decoding) mass spectra
# using multiprocessing.
_UNPACK_BATCH_SIZE = 64


class UReader():
    """
    Minimal object for reading data from UIMF files
    """

    def __init__(self, path, workers):
        """
        _UReader initialized with path to UIMF file. Should be closed when done

        TODO (Dylan Ross): accomodate manually modifying the mass calibration parameters

        Parameters
        ----------
        path : ``str``
            path to UIMF file
        workers : ``int``
            number of worker processes for deconding and decompressing spectra
        """
        self._con = connect(path)
        self._cur = self._con.cursor()
        # store number of m/z bins
        self._n_bins = self._get_n_bins()
        # store number of frames
        self._n_frames = self._get_n_frames()
        # store the max number of scans per frame (Prescan_TOFPulses)
        self._max_frame_n_scans = self._get_max_frame_n_scans()
        # store the average duration of all IM frames (in milliseconds)
        self._avg_frame_duration = self._get_avg_frame_duration()
        # store array of all m/zs
        self._all_mzs = self.mzbin_to_mz(np.arange(self._n_bins))
        # set up multiprocessing queues for decoding and decompressing spectra
        self._q_in = JoinableQueue()
        # NOTE (Dylan Ross): The real limiting factor as discussed in the note above
        #                    has more to do with the size of items in the output queue
        #                    rather than the input queue since the spectra get decompressed
        #                    and decoded into much larger datastructures. In this case it 
        #                    help to use separate output queues for each of the worker
        #                    threads so that added size gets spread over n_threads chunks
        #                    increasing the batch size that can be successfully processed 
        #                    without running into problems
        self._q_outs = [
            JoinableQueue() for n in range(workers)
        ]
        # set up worker processes for decoding and decompressing spectra
        self._workers = [
            Process(target=self._worker, name=f"worker_{n}", 
                    args=(self._q_in, self._q_outs[n]), 
                    daemon=True) 
                for n in range(workers)
        ]
        # start workers for decoding and decompressing spectra
        for worker in self._workers:
            worker.start()
        # init shared memory for accumulated spectrum
        # https://docs.python.org/3/library/multiprocessing.shared_memory.html#:~:text=the%20shared%20memory-,The%20following%20example,-demonstrates%20a%20practical
        #accum_spec = np.zeros(self._n_bins)
        #self.shm = shared_memory.SharedMemory(name="accum_spec_mem", create=True, size=accum_spec.nbytes)
        

    @staticmethod
    def _worker(q_in, q_out):
        """
        process worker function for decoding and decompressing spectra
        """
        n_processed = 0
        pid = os.getpid()
        while True:
            # block, wait for a new job from the input queue
            #print(f"(pid: {pid}) waiting for input...")
            tag, blob = q_in.get()        
            #print(f"(pid: {pid}) got input (tag: {tag})")
            spectrum = UReader._decode(UReader._decompress(blob))
            #print(f"(pid: {pid}) unpacked spectrum, adding to output queue (tag: {tag})")
            q_out.put((tag, spectrum))
            # signal to input queue that we finished our job
            q_in.task_done()
            n_processed += 1
            #print(f"(pid: {pid}) processing complete (tag: {tag})")
            #print(f"(pid: {pid}) processed {n_processed} spectra")

    def close(self):
        """ clean up tasks """
        # close DB connection
        self._con.close()
        # terminate all of the worker processes
        for worker in self._workers:
            worker.terminate()

    def _get_n_bins(self):
        """ fetches NBins value from Global_Params table """
        return int(self._cur.execute("SELECT ParamValue FROM Global_Params WHERE ParamName='Bins'").fetchone()[0])
    
    def _get_n_frames(self):
        """ fetches N_frames value from Global_Parameters table """
        return int(self._cur.execute('SELECT NumFrames FROM Global_Parameters').fetchone()[0])

    def _get_max_frame_n_scans(self):
        """ fetch the max n_scans for all frames (Prescan_TOFPulses from Global_Parameters) """
        return int(self._cur.execute('SELECT Prescan_TOFPulses FROM Global_Parameters').fetchone()[0])

    def _get_avg_frame_duration(self):
        """ computes the average frame duration across all IM frames in milliseconds """
        durations = [_[0] for _ in self._cur.execute('SELECT Duration FROM Frame_Parameters').fetchall()]
        return np.mean(durations) * 1000  # multiply by 1000 for milliseconds 

    def scan_to_ms(self, scan):
        """ 
        converts scan index into scan time in milliseconds, 

        Parameters
        ----------
        scan : ``int``
            scan index, operation is broadcastable to numpy.ndarray(int) for multiple scans if desired

        Returns
        -------
        scan_time : ``float``
            scan time in milliseconds
        """
        return  self._avg_frame_duration * scan / self._max_frame_n_scans

    @staticmethod
    def _decompress(blob):
        """
        Use LZF to decompress the encoded intensity array

        Parameters
        ----------
        blob : ``bytes``
            compressed and encoded intensity array, directly from SQLite3 query

        Returns
        -------
        enc_array : ``numpy.array(int32)``
            encoded intensity array (as 32 bit ints)
        """
        l = len(blob)
        if l == 0:
            return np.array([], dtype=np.int32)
        else:
            return np.frombuffer(lzf.decompress(blob, 100 * l), dtype=np.int32)

    @staticmethod
    def _decode(enc_array):
        """
        Decodes the encoded intensity array (output from _decompress), produces array of mzbins and intensities

        TODO (Dylan Ross): What type of encoding is being used here? -> RLZE!

        Parameters
        ----------
        enc_array : ``numpy.array(int32)``
            encoded intensity array (as 32 bit ints)

        Returns
        -------
        mzb_i : ``list(tuple(int, int))``
            pairs of m/z bin and intensity values
        """
        mzb_i = []
        #prev_value = 0
        bin_idx = 0
        for value in enc_array:
            if value < 0:
                # negative values increment bin index
                bin_idx -= value
            #elif value == 0 and prev_value < -1e10:
            #    # ?
            #    pass
            else:
                mzb_i.append((bin_idx, value))
                bin_idx += 1
                #prev_value = value
        return mzb_i
    
    @staticmethod
    def _decode_for_atd(enc_array, mzbin_min, mzbin_max):
        i_sum = 0
        bin_idx = 0
        for value in enc_array:
            if value < 0:
                # negative values increment bin index
                bin_idx -= value
            else:
                if bin_idx >= mzbin_min:
                    if bin_idx <= mzbin_max:
                        i_sum += value
                        bin_idx += 1
                    else:
                        break
        return i_sum
    
    def _unpack_batch(self, tags_blobs):
        """
        takes a list of tuples (tag, blob) with compressed encoded spectra, uses
        multiprocessing to decompress and decode them, returning a list of tuples
        (tag, intensities) with decompressed and decoded spectra  
        """
        # do not accept too many inputs, the errors are difficult to debug
        assert(len(tags_blobs) <= _UNPACK_BATCH_SIZE)
        for tag, blob in tags_blobs:
            self._q_in.put((tag, blob))
        # NOTE (Dylan Ross): this will stop blocking after all items that were retrieved have had
        #                    task_done() called in the consumer, marking them as complete. The 
        #                    task_done() method is only used in the consumer after the output has 
        #                    been put into the output queue so we can trust that all processing is
        #                    complete once the input queue finally stops blocking
        self._q_in.join()
        # fetch the spectra from the output queue
        tags_intensities = []
        # iterate through each thread output queue
        for q_out in self._q_outs:
            while True:  # emulate a do-while loop
                try:
                    tags_intensities.append(q_out.get_nowait())
                except queue.Empty:
                    break
        return tags_intensities

    def _accum_spectra(self, qry):
        """
        backend method for accumulating spectra, takes a query then does the extraction and processing 
        and returns the intensities array. Enables using slightly different queries without repeating the 
        logic for the actual data extraction and processing.

        Parameters
        ----------
        qry : ``str``
            query for selecting spectra (by frame and/or scan)
        
        Returns
        -------
        mzs : ``numpy.array(float)``
            m/z values for summed mass spectrum
        intensities : ``numpy.array(int)``
            summed mass spectrum, indices are m/z bins
        """
        res = self._cur.execute(qry).fetchall()
        intensities = np.zeros(self._n_bins)
        batch_n = 0
        batches = len(res) // _UNPACK_BATCH_SIZE + 1
        for i in range(0, len(res), _UNPACK_BATCH_SIZE):
            batch_n += 1
            print(f"\rdecoding spectra batch: {batch_n:6d}/{batches:<6d} ({100. * batch_n / batches:6.2f} %)", end=" ")
            batch = res[i:i + _UNPACK_BATCH_SIZE]
            for tag, spectrum in self._unpack_batch(batch):
                for mzb, i in spectrum:
                    intensities[mzb] += i
        print()
        return self._all_mzs, intensities

    def accum_spectra(self, frames, scans):
        """ 
        sum together spectra from multiple frames/scans, returns array of summed intensities by m/z bins
        works best for targeting a small number of frames/scans

        Parameters
        ----------
        frames : ``list(int)``
        scans : ``list(int)``
            lists of frames/scans to accumulate spectra from

        Returns
        -------
        intensities : ``numpy.array(int)``
            summed mass spectrum, indices are m/z bins
        """
        f = ','.join([str(_) for _ in frames])
        s = ','.join([str(_) for _ in scans])
        qry = "SELECT ScanNum, Intensities FROM Frame_Scans WHERE FrameNum IN ({}) AND ScanNum IN ({})"
        return self._accum_spectra(qry.format(f, s))

    def accum_spectra_sparsescans(self, frames, scan_mod):
        """ 
        sum together spectra from specified frames and sparsely sampled scans, returns array of summed intensities by m/z bins
        include only 1 in {scan_mod} scans
        works best for targeting a small number of frames

        Parameters
        ----------
        frames : ``list(int)``
            list of frames to accumulate spectra from
        scan_mod : ``int``
            determines how densely the arrival time dimension is sampled, 1 in {scan_mod} scans are included

        Returns
        -------
        intensities : ``numpy.array(int)``
            summed mass spectrum, indices are m/z bins
        """
        f = ','.join([str(_) for _ in frames])
        qry = "SELECT ScanNum, Intensities FROM Frame_Scans WHERE FrameNum IN ({}) AND ScanNum%{}==0".format(f, scan_mod)
        return self._accum_spectra(qry)

    def accum_spectra_allframes(self, 
                                scans=None, scan_min=None, scan_max=None, skip_frame_1=False):
        """ 
        sum together spectra from all frames and specified scans, returns array of summed intensities by m/z bins
        scans can be specified as a list or as min/max values for a range

        Parameters
        ----------
        scans : ``list(int)``, optional
            list of scans to accumulate spectra from
        scan_min, scan_max : ``int``, optional
            alternative to scans, specify min/max scan values to include a range of scans
            ignored if scans is provided or if only min or max is provided
        skip_frame_1 : ``bool``, default=False
            skip the first frame when accumulating spectra from multiple frames

        Returns
        -------
        intensities : ``numpy.array(int)``
            summed mass spectrum, indices are m/z bins
        """
        qry = "SELECT ScanNum, Intensities FROM Frame_Scans"
        if scans is not None:
            s = ','.join([str(_) for _ in scans])
            qry += " WHERE ScanNum IN ({})".format(s)
        elif scan_min is not None and scan_max is not None:
            qry += " WHERE ScanNum>={} AND ScanNum<={}".format(scan_min, scan_max)
        if skip_frame_1:
            where = "WHERE " if "WHERE" not in qry else "AND "
            qry += " " + where + "FrameNum>1"
        return self._accum_spectra(qry)
    
    def extract_atd_allframes(self,
                              target_mz, mz_ppm,
                              skip_frame_1=False):
        """
        Extract an arrival time distribution for a target m/z +/- ppm
        """
        # compute the m/z tolerance from ppm
        tol = target_mz * mz_ppm / 1e6
        # figure out the min/max mz bins to include
        mzbin_min = self.mz_to_mzbin(target_mz - tol)
        mzbin_max = self.mz_to_mzbin(target_mz + tol)
        # set up the query to get the data
        qry = "SELECT ScanNum, Intensities FROM Frame_Scans" 
        if skip_frame_1:
            qry += " WHERE FrameNum>1"
        # create the ATD arrays to accumulate into
        #atd_at = self.scan_to_ms(np.arange(self._max_frame_n_scans))
        atd_i = np.zeros(self._max_frame_n_scans)
        for scan, blob in self._cur.execute(qry):
            atd_i[scan] += UReader._decode_for_atd(UReader._decompress(blob), 
                                                    mzbin_min, mzbin_max)
        return self.scan_to_ms(np.arange(self._max_frame_n_scans)), atd_i        
    
    def _accum_sparse_data(self, 
                           qry: str, 
                           i_min: int, 
                           include_scans: Optional[Set[int]] = None) :
        """
        accumulate sparse data retaining m/z and arrival time scans

        Parameters
        ----------
        qry : ``str``
            query for selecting spectra (by frame and/or scan)
        i_min : ``int``
            minimum intensity to include points
        include_scans : ``set(int)``, optional
            if provided, limit the data that is extracted and decoded to 
            only include scan numbers that are present in the set
        
        Returns
        -------
        scans : ``numpy.ndarray(int)``
            1D array of scans
        mzs : ``numpy.ndarray(float)``
            1D array of m/zs
        intensites : ``numpy.ndarray(int)``
            1D array of intensities 
        """
        res = self._cur.execute(qry).fetchall()
        if include_scans is not None:
            res = [_ for _ in res if _[0] in include_scans]
        n = 0
        batch_n = 0
        batches = len(res) // _UNPACK_BATCH_SIZE + 1
        scans, mzs, intensities = [], [], []
        for i in range(0, len(res), _UNPACK_BATCH_SIZE):
            batch_n += 1
            print(f"\rdecoding spectra batch: {batch_n:6d}/{batches:<6d} ({100. * batch_n / batches:6.2f} %)", end=" ")
            batch = res[i:i + _UNPACK_BATCH_SIZE]
            for scan, spectrum in self._unpack_batch(batch): 
                for mzb, i in spectrum:
                    if i >= i_min:
                        scans.append(scan)
                        mzs.append(self._all_mzs[mzb])
                        intensities.append(i)
                        # keep track of the number of points included
                        n += 1
        print()
        # construct the return arrays
        # NOTE: https://stackoverflow.com/questions/30012362/faster-way-to-convert-list-of-objects-to-numpy-array
        sparse_data = [
            np.fromiter(scans, np.int32, n), 
            np.fromiter(mzs, np.float64, n),
            np.fromiter(intensities, np.int32, n)
        ]
        return sparse_data

    def extract_survey_data_sparse(self, frame_mod, scan_mod, i_min):
        """
        extracts and decodes spectra
        from 1 out of every {frame_mod} frames and 
        1 out of every {scan_mod} scans selected
        ! skips first frame !
        maintains spectra as sparse representations 
        
        Used for targeted extraction, this produces survey data to check for targets 
        and guide full data extraction for targets in subsequent steps.
        
        Parameters
        ----------
        frame_mod : ``int``
            determines how densely the frames are sampled in survey data
            1 in {frame_mod} frames are selected
        scan_mod : ``int``
            determines how densely the arrival time (scans) is sampled in survey data
            1 in {scan_mod} scans are selected
        i_min : ``float``
            minimum intensity cutoff for individual spectrum points
            
        Returns
        -------
        survey_scans : ``numpy.ndarray(int)``
            1D array of scans
        survey_mzs : ``numpy.ndarray(float)``
            1D array of m/zs
        survey_intensites : ``numpy.ndarray(int)``
            1D array of intensities 
        """
        qry = ('SELECT ScanNum SORTED, Intensities FROM Frame_Scans '
               'WHERE FrameNum%{}==0 AND FrameNum>1 AND ScanNum>0 AND ScanNum%{}==0')
        return self._accum_sparse_data(qry.format(frame_mod, scan_mod), i_min)
    
    def extract_full_data_sparse(self, 
                                 i_min: float,
                                 include_scans: Optional[Set[int]] = None) :
        """
        extracts and decodes spectra from selected scans and all frames
        ! skips first frame !
        maintains spectra as sparse representations 
        
        Used for targeted extraction, this produces survey data to check for targets 
        and guide full data extraction for targets in subsequent steps.
        
        Parameters
        ----------
        i_min : ``float``
            minimum intensity cutoff for individual spectrum points
        include_scans : ``set(int)``, optional
            if provided, limit the data that is extracted and decoded to 
            only include scan numbers that are present in the set
        
        Returns
        -------
        scans : ``numpy.ndarray(int)``
            1D array of scans
        mzs : ``numpy.ndarray(float)``
            1D array of m/zs
        intensites : ``numpy.ndarray(int)``
            1D array of intensities 
        """
        qry = 'SELECT ScanNum SORTED, Intensities FROM Frame_Scans WHERE FrameNum > 1'
        return self._accum_sparse_data(qry, i_min, include_scans=include_scans)
    
    def extract_full_data_sparse_scan_range(self, 
                                            i_min: float,
                                            scan_min: int, 
                                            scan_max: int) :
        """
        extracts and decodes spectra from selected scans and all frames
        ! skips first frame !
        maintains spectra as sparse representations 
        
        Used for targeted extraction, this produces survey data to check for targets 
        and guide full data extraction for targets in subsequent steps.
        
        Parameters
        ----------
        i_min : ``float``
            minimum intensity cutoff for individual spectrum points
        scan_min : ``int``
        scan_max : ``int``
            restrict included scans to a defined range
        
        Returns
        -------
        scans : ``numpy.ndarray(int)``
            1D array of scans
        mzs : ``numpy.ndarray(float)``
            1D array of m/zs
        intensites : ``numpy.ndarray(int)``
            1D array of intensities 
        """
        qry = """--sqlite3
        SELECT 
            ScanNum SORTED, 
            Intensities 
        FROM 
            Frame_Scans 
        WHERE 
            FrameNum > 1
            AND ScanNum >= {smin}
            AND ScanNum <= {smax}
        ;""".format(smin=scan_min, smax=scan_max)
        return self._accum_sparse_data(qry, i_min)
    
    def _get_mz_cal_params(self, frames):
        """ get m/z calibration parameters for a specified range of frames, average them together """
        qry = "SELECT ParamValue FROM V_Frame_Params WHERE ParamName='{}'"
        slope_qry = qry.format('CalibrationSlope')
        intercept_qry = qry.format('CalibrationIntercept')
        if frames is not None:
            f = ','.join([str(_) for _ in frames])
            slope_qry += " AND FrameNum IN ({})".format(f)
            intercept_qry += " AND FrameNum IN ({})".format(f)
        slopes = [float(_[0]) for _ in self._cur.execute(slope_qry).fetchall()]
        intercepts = [float(_[0]) for _ in self._cur.execute(intercept_qry).fetchall()]
        return np.mean(slopes), np.mean(intercepts)
    
    def _get_mz_cal_params_all_frames(self):
        """ get m/z calibration parameters for all frames, average them together """
        # set frames to None to include all frames
        return self._get_mz_cal_params(None)

    def _get_bin_width(self):
        """ fetches BinWidth value from GlobalParams table """
        qry = "SELECT ParamValue FROM Global_Params WHERE ParamName='BinWidth'"
        return float(self._cur.execute(qry).fetchone()[0])
    
    def mzbin_to_mz(self, mzbin):
        """
        convert mzbin to m/z using TOF calibration parameters

        Parameters
        ----------
        mzbin : ``int``
            mzbin

        Returns
        -------
        mz : ``float``
            m/z
        """
        bin_width = self._get_bin_width()
        cal_slope, cal_intercept = self._get_mz_cal_params_all_frames()
        return (cal_slope / 1e4 * (mzbin * bin_width * 10. - cal_intercept * 1e4))**2.
    
    def mz_to_mzbin(self, mz):
        # this is really terrible, there is a much better way to do this to be sure
        return np.argmin(np.abs(self._all_mzs - mz))


class UReaderSetMzParams(UReader):
    """  """

    def __init__(self, path, workers, mz_cal_slope, mz_cal_intercept):
        self.mz_cal_slope, self.mz_cal_intercept = mz_cal_slope, mz_cal_intercept
        super().__init__(path, workers)
        

    def _get_mz_cal_params_all_frames(self):
        return self.mz_cal_slope, self.mz_cal_intercept
