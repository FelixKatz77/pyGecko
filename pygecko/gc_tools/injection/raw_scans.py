from typing import Iterator

import numpy as np
import pandas as pd


class Raw_Scans:

    '''
    Class to hold the centroids of an MS injection exactly as they were read.

    The layout is the flat one ANDI/AIA uses: one m/z and one intensity array for the whole run,
    plus the offset at which each scan starts. m/z is kept unrounded and the intensity array keeps
    the dtype the source delivered, so an export can write back what was read. The nominal-mass
    matrix the analysis works on is derived from here, in one place, by to_nominal_matrix.

    Attributes:
        retention_times (np.ndarray): Start time of each scan in milliseconds.
        scan_index (np.ndarray): Offset of each scan's first centroid in mz and intensity.
        mz (np.ndarray): m/z of every centroid of the run, unrounded, in source order.
        intensity (np.ndarray): Intensity of every centroid of the run, in the source dtype.
    '''

    retention_times: np.ndarray
    scan_index: np.ndarray
    mz: np.ndarray
    intensity: np.ndarray

    __slots__ = 'retention_times', 'scan_index', 'mz', 'intensity'

    def __init__(self, retention_times: np.ndarray, scan_index: np.ndarray, mz: np.ndarray,
                 intensity: np.ndarray):
        self.retention_times = np.asarray(retention_times, dtype=np.float64)
        self.scan_index = np.asarray(scan_index, dtype=np.int64)
        self.mz = np.asarray(mz, dtype=np.float64)
        self.intensity = np.asarray(intensity)

    def spectra(self) -> Iterator[tuple[float, np.ndarray, np.ndarray]]:

        '''
        Yields the retention time in milliseconds, the m/z array and the intensity array of every
        scan in source order.
        '''

        bounds = np.append(self.scan_index, len(self.mz))
        for retention_time, start, stop in zip(self.retention_times, bounds[:-1], bounds[1:]):
            yield retention_time, self.mz[start:stop], self.intensity[start:stop]

    def to_nominal_matrix(self) -> pd.DataFrame:

        '''
        Returns the scans binned to nominal mass: one row per scan indexed by retention time in
        milliseconds, one sorted integer column per nominal m/z, the summed intensity where a
        scan has several centroids rounding to the same nominal mass, and 0 where it has none.

        Summing conserves the ion count, so the TIC derived from the matrix equals the sum of the
        centroids as acquired.
        '''

        nominal_masses = np.round(self.mz).astype(int)
        scan_lengths = np.diff(np.append(self.scan_index, len(self.mz)))
        rows = np.repeat(np.arange(len(self.retention_times)), scan_lengths)
        columns, column_indices = np.unique(nominal_masses, return_inverse=True)
        matrix = np.zeros((len(self.retention_times), len(columns)))
        np.add.at(matrix, (rows, column_indices), self.intensity.astype(np.float64))
        return pd.DataFrame(matrix, columns=columns,
                            index=pd.Index(self.retention_times, name='retention_time'))
