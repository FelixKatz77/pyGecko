from typing import Sequence
import numpy as np

class Utilities:

    @staticmethod
    def convert_time_to_scan(time_value: float|Sequence[float], scan_rate: float) -> list[int]:

        '''
        Takes in single time or sequence of two times in minutes, returns the closest scan index or a sequence of two
        indices.
        '''

        if not type(time_value) in [np.ndarray, list, tuple]: time_value = [time_value]
        indices = []
        for value in time_value:
            value = int(np.round(value / scan_rate))
            indices.append(value)
        if len(indices) == 1: indices = indices[0]
        return indices

    @staticmethod
    def check_interval(value:float, midpoint:float, tolerance:float) -> bool:

        '''
        Takes in a value, midpoint, and tolerance and returns True if the value is within the interval defined by the
        midpoint and tolerance.
        '''

        if (midpoint-tolerance) < value < (midpoint+tolerance):
            return True
        else:
            return False

    @staticmethod
    def find_empty_ranges(chromatogram: np.ndarray, *, threshold: float = 0.0,
                          min_duration: float = 0.0) -> list[tuple[float, float]]:
        '''
        Find contiguous time ranges where the signal is "empty" (<= threshold or NaN).

        Args:
            time_axis (np.ndarray): 1D array of times (same length as signal), monotonically increasing.
            signal (np.ndarray): 1D array of signal values.
            threshold (float): Values <= threshold are considered empty. Default 0.0.
            min_duration (float): Minimum duration (in time units of time_axis) to report a gap. Default 0.0.

        Returns:
            list[tuple[float, float]]: List of (start_time, end_time) for empty ranges.
        '''

        time_axis, signal = chromatogram

        if time_axis.ndim != 1 or signal.ndim != 1:
            raise ValueError("time_axis and signal must be 1D arrays")
        if len(time_axis) != len(signal):
            raise ValueError("time_axis and signal must have the same length")
        if len(time_axis) == 0:
            return []

        # Empty where NaN or <= threshold
        empty_mask = np.isnan(signal) | np.any(signal == 0)

        ranges: list[tuple[float, float]] = []
        if not np.any(empty_mask):
            return ranges

        # Find transitions in the mask
        # We pad with False at both ends to catch edges cleanly
        padded = np.pad(empty_mask.astype(np.int8), (1, 1), constant_values=0)
        diff = np.diff(padded)

        # Starts where diff == 1, ends where diff == -1 (end index is exclusive)
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]

        for s_idx, e_idx in zip(starts, ends):
            # Convert index range [s_idx, e_idx) to time range [start_time, end_time]
            # Use inclusive end index e_idx - 1 to fetch end time present in data
            start_time = float(time_axis[s_idx])
            end_time = float(time_axis[e_idx - 1])
            if end_time < start_time:
                continue
            duration = end_time - start_time
            if duration + 0.0 >= min_duration:
                ranges.append((start_time, end_time))

        return ranges

