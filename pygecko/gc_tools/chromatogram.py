from typing import Literal

import numpy as np

from pygecko.gc_tools.utilities import Utilities


class Chromatogram:

    '''
    Detector signal of one injection on a minutes time axis.

    Peaks are not held here: they stay on Injection.peaks, where matching, flagging and RI
    assignment operate and where MS peaks carry the spectra read from the injection's scans.

    Attributes:
        time (np.ndarray): Retention time of each scan in minutes, strictly increasing.
        intensity (np.ndarray): Raw detector signal, one value per scan.
        processed (np.ndarray|None): Smoothed, baseline-corrected signal on the same axis as time; None until
        a baseline correction has run.
        kind (str): 'FID' for a flame ionization trace, 'TIC' for an MS total ion current.
    '''

    time: np.ndarray
    intensity: np.ndarray
    processed: np.ndarray|None
    kind: Literal['FID', 'TIC']

    __slots__ = 'time', 'intensity', 'processed', 'kind'

    def __init__(self, time: np.ndarray, intensity: np.ndarray, kind: Literal['FID', 'TIC'],
                 processed: np.ndarray|None = None):
        time = np.asarray(time, dtype=np.float64)
        intensity = np.asarray(intensity)
        if time.ndim != 1 or time.shape != intensity.shape:
            raise ValueError(f'time and intensity must be 1D arrays of equal length, got {time.shape} and '
                             f'{intensity.shape}.')
        if np.any(np.diff(time) <= 0):
            raise ValueError('time must be strictly increasing.')
        self.time = time
        self.intensity = intensity
        self.processed = processed
        self.kind = kind

    def __getitem__(self, index: slice) -> 'Chromatogram':

        '''
        Returns the scans in index as a new Chromatogram, with the processed signal sliced alongside.
        '''

        processed = None if self.processed is None else self.processed[index]
        return Chromatogram(self.time[index], self.intensity[index], self.kind, processed)

    @property
    def scan_rate(self) -> float:

        '''
        Returns the sampling interval in minutes, taken between the second and third scan as
        Analysis_Settings did before the chromatogram carried it.
        '''

        return float(self.time[2] - self.time[1])

    @property
    def run_time(self) -> tuple[float, float]:

        '''
        Returns the first and last retention time of the chromatogram in minutes.
        '''

        return float(self.time[0]), float(self.time[-1])

    def empty_ranges(self) -> list[tuple[float, float]]:

        '''
        Returns the (start, end) retention times in minutes of every stretch without signal.

        A real MS TIC carries dozens of short dropouts per injection, which is why they are reported
        on request here rather than announced when an injection is built.
        '''

        return Utilities.find_empty_ranges((self.time, self.intensity))
