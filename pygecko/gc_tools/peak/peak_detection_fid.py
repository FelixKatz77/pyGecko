import numpy as np
from numpy import ndarray
from pybaselines import Baseline
from scipy.integrate import simpson
from scipy.signal import find_peaks, savgol_filter, argrelmin
from scipy.ndimage import gaussian_filter1d
from statsmodels.stats.stattools import durbin_watson
from copy import copy
from pygecko.gc_tools.analysis.analysis_settings import Analysis_Settings
from pygecko.gc_tools.chromatogram import Chromatogram
from pygecko.gc_tools.peak.fid_peak import FID_Peak
from pygecko.gc_tools.utilities import Utilities


class Peak_Detection_FID:
    '''
    A class wrapping functions to detect peaks in FID chromatograms.
    '''

    @staticmethod
    def baseline_correction(chromatogram: Chromatogram, analysis_settings: Analysis_Settings) -> np.ndarray:

        '''
        Returns the baseline corrected signal of a chromatogram.

        The whole chromatogram is corrected, whatever time_range is set to: the window only selects
        which peaks pick_peaks keeps. Cropping here made the processed signal, the savgol window
        auto-tuning and the SNIP edge effects depend on the window, and a later wider window could
        not recover signal the crop had cut away.

        Args:
            chromatogram (Chromatogram): Chromatogram to correct.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the baseline correction.

        Returns:
            np.ndarray: Smoothed, baseline corrected intensity on chromatogram.time.
        '''

        y_smooth = Peak_Detection_FID.__savgol(chromatogram.intensity, analysis_settings)
        y_corr, baseline = Peak_Detection_FID.__baseline_filter(chromatogram.time, y_smooth, analysis_settings)
        return y_corr

    @staticmethod
    def pick_peaks(chromatogram: Chromatogram, analysis_settings: Analysis_Settings) -> dict[float:FID_Peak]:

        '''
        Returns a dictionary of FID peaks detected in the processed signal of a chromatogram.

        Args:
            chromatogram (Chromatogram): Chromatogram with a processed signal to detect peaks in.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            dict[float:FID_Peak]: Dictionary of FID peaks.
        '''

        peak_rts, peak_widths, peak_heights, peak_boarders, peak_areas, flag_peaks_list = Peak_Detection_FID.__detect_peaks(
            chromatogram.time, chromatogram.processed, chromatogram.scan_rate, analysis_settings)
        peaks = Peak_Detection_FID.__initialize_peaks(peak_rts, peak_heights, peak_widths, peak_boarders, peak_areas, flag_peaks_list)
        return peaks



    @staticmethod
    def peak_area(time: np.ndarray, signal: np.ndarray, start: int, end: int) -> float:

        '''
        Returns the area of a signal between two scan indices in intensity x minutes.

        Integrating over the time axis rather than over scans makes the area independent of the
        sampling rate. pick_peaks and FID_Injection.integrate both call this, so they agree.

        Args:
            time (np.ndarray): Retention times of the signal in minutes.
            signal (np.ndarray): Signal to integrate.
            start (int): Index of the first scan of the peak.
            end (int): Index one past the last scan of the peak.

        Returns:
            float: Area of the peak.
        '''

        return float(simpson(signal[start:end], x=time[start:end]))

    @staticmethod
    def __detect_peaks(time: np.ndarray, signal: np.ndarray, scan_rate: float, analysis_settings: Analysis_Settings):

        '''
        Returns the peak retention times, widths, heights, boarders and areas of a chromatogram.

        Peaks are detected over the whole signal; a time_range then keeps those whose apex lies
        inside it. Their boarders may reach past the window, so an integration is never cut in half,
        and a peak kept by a window is identical to the same peak picked over the full run.

        Args:
            time (np.ndarray): Retention times of the signal in minutes.
            signal (np.ndarray): Baseline corrected signal to detect peaks in.
            scan_rate (float): Sampling interval of the signal in minutes.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[int]]: Peak retention times,
            widths, heights, boarders, areas and the indices of the peaks to flag for boarder overlap.
        '''

        prominence = analysis_settings.pop('prominence_fid', np.mean(signal))
        width = analysis_settings.pop('width', 0)
        #TODO: Implement S/N.
        height = analysis_settings.pop('height', 0)

        peak_indices, peak_properties = find_peaks(signal, prominence=prominence, width=width,
                                                   height=height)
        peak_widths, peak_heights = peak_properties['widths'], peak_properties['peak_heights']
        peak_boarders = Peak_Detection_FID.__detect_borders(signal, peak_indices, peak_widths, scan_rate,
                                                             analysis_settings)
        peak_boarders, flag_peaks_list = Peak_Detection_FID.__resolve_boarder_overlap(peak_boarders, peak_indices, signal)
        peak_areas = Peak_Detection_FID.__calculate_areas(time, signal, peak_boarders)
        peak_widths = np.array([((boarder[1] - boarder[0])*scan_rate) for boarder in peak_boarders])
        # Boarders and indices index the signal, so their retention times are read straight off its
        # own time axis rather than reconstructed from scan_rate: index * scan_rate + t0 lands about
        # 1e-15 from the stored axis value, enough for searchsorted in FID_Injection.integrate to pick
        # the neighbouring scan. __find_right_boarder may return one past the last scan, a valid slice
        # bound but not a valid time, so the right boarder is clamped to the last scan.
        peak_boarders = time[np.clip(peak_boarders, 0, len(time) - 1)]
        peak_rts = time[peak_indices]

        time_range = analysis_settings.pop('time_range', None)
        if time_range:
            kept = np.flatnonzero((peak_rts >= time_range[0]) & (peak_rts <= time_range[1]))
        else:
            kept = np.arange(len(peak_rts))
        flag_peaks_list = [j for j, i in enumerate(kept) if i in flag_peaks_list]
        return (peak_rts[kept], peak_widths[kept], peak_heights[kept], peak_boarders[kept],
                np.asarray(peak_areas)[kept], flag_peaks_list)

    @staticmethod
    def __baseline_filter(time: np.ndarray, signal: np.ndarray, analysis_settings: Analysis_Settings) -> tuple[np.ndarray, np.ndarray]:

        '''
        Returns a baseline corrected signal and the baseline.

        Args:
            time (np.ndarray): Retention times of the signal in minutes.
            signal (np.ndarray): Signal to correct.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the baseline correction.

        Returns:
            tuple[np.ndarray, np.ndarray]: Baseline corrected signal and baseline.
        '''

        max_half_window = analysis_settings.pop('max_half_window', 200)
        baseline_fitter = Baseline(x_data=time)
        baseline = baseline_fitter.snip(signal, max_half_window=max_half_window)[0]
        y_corr = signal - baseline
        return y_corr, baseline

    @staticmethod
    def __savgol(signal: np.ndarray, analysis_settings: Analysis_Settings) -> np.ndarray:

        '''
        Returns a Savitzky-Golay filtered signal.

        Args:
            signal (np.ndarray): Signal to filter.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the Savitzky-Golay filter.

        Returns:
            np.ndarray: Savitzky-Golay filtered signal.
        '''

        savgol_window = analysis_settings.pop('savgol_window', Peak_Detection_FID.__optimize_savgol_window(signal))
        return savgol_filter(signal, savgol_window, 2)

    @staticmethod
    def __optimize_savgol_window(signal: np.ndarray) -> int:

        '''
        Returns al.

        Args:
            signal: Signal to optimize the Savitzky-Golay window size for.

        Returns:
            int: Optimal Savitzky-Golay window size.
        '''

        windows = [5, 7, 9, 11, 21, 31]
        best_dw = 0
        for window in windows:
            _y = savgol_filter(signal, window, 2)
            resids = signal - _y
            dw = durbin_watson(resids)
            if abs(2 - dw) < abs(2 - best_dw):
                best_dw = dw
                best_window = window
        return best_window

    @staticmethod
    def __detect_borders(signal: np.ndarray, peak_indices: np.ndarray, peak_widths: np.ndarray, scan_rate: float,
                          analysis_settings: Analysis_Settings) -> np.ndarray:

        '''
        Returns the boarders of the peaks in a chromatogram.

        Args:
            signal (np.ndarray): Signal to detect the boarders for.
            peak_indices (np.ndarray): Indices of the signal at which peaks are located.
            peak_widths (np.ndarray): Widths of the peaks.
            scan_rate (float): Sampling interval of the signal in minutes.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the boarder detection.

        Returns:
            np.ndarray: Boarders of the peaks.
        '''

        first_diff = np.diff(signal) / scan_rate
        boarder_threshold = analysis_settings.pop('boarder_threshold', abs(np.mean(first_diff))*0.5)
        boarder_window = analysis_settings.pop('boarder_window', 100)
        # Prefix sums let the boarder searches read a window mean as one subtraction. The searches
        # step one scan at a time and averaged the window afresh at every step, which was three
        # quarters of pick_peaks on a 30k-scan chromatogram.
        diff_cumsum = np.concatenate(([0.0], np.cumsum(first_diff)))
        boarders = np.empty((len(peak_indices), 2), int)

        for i, index in enumerate(peak_indices):
            left = Peak_Detection_FID.__find_left_boarder(index, peak_widths[i], diff_cumsum, boarder_threshold,
                                                          boarder_window)
            right = Peak_Detection_FID.__find_right_boarder(signal, index, peak_widths[i], diff_cumsum,
                                                            boarder_threshold, boarder_window)
            boarders[i] = [left, right]
        return boarders

    @staticmethod
    def __find_left_boarder(index:int, peak_width:int, diff_cumsum:np.ndarray, threshold:float, window:int) -> int:

        '''
        Returns the left boarder of a peak in a chromatogram.

        Args:
            index (int): The index of the peak to find the left boarder for.
            peak_width (int): The width of the peak.
            diff_cumsum (np.ndarray): Prefix sums of the first derivative of the chromatogram, with a
            leading 0, so that diff_cumsum[b] - diff_cumsum[a] is the sum of first_diff[a:b].
            threshold (float): The threshold for the first derivative.
            window (int): The window size for the boarder detection.

        Returns:
            int: The left boarder of the peak.
        '''

        _l = int(round(index - peak_width / 2, 0))
        if _l <= 0:
            return 0
        left = max(0, _l - window)

        while (diff_cumsum[_l] - diff_cumsum[int(left)]) / (_l - int(left)) > threshold:
            _l -= 1
            if _l <= 0:
                break
            left = max(0, _l - window)
        left = max(0, _l - window / 4)
        return left

    @staticmethod
    def __find_right_boarder(chromatogram: np.ndarray, index: int, peak_width: int, diff_cumsum: np.ndarray,
                             threshold: float, window: int) -> int:

        '''
        Returns the right boarder of a peak in a chromatogram.

        Args:
            chromatogram (np.ndarray): The intensity values of the chromatogram.
            index (int): The index of the peak to find the right boarder for.
            peak_width (int): The width of the peak.
            diff_cumsum (np.ndarray): Prefix sums of the first derivative of the chromatogram, with a
            leading 0, so that diff_cumsum[b] - diff_cumsum[a] is the sum of first_diff[a:b].
            threshold (float): The threshold for the first derivative.
            window (int): The window size for the boarder detection.

        Returns:
            int: The right boarder of the peak.
        '''

        diff_length = diff_cumsum.shape[0] - 1
        _r = int(round(index + peak_width / 2, 0))
        if _r >= diff_length:
            return chromatogram.shape[0]
        right = min(diff_length, _r + window)

        while abs((diff_cumsum[int(right)] - diff_cumsum[_r]) / (int(right) - _r)) > threshold:
            _r += 1
            if _r >= diff_length:
                break
            right = min(diff_length, _r + window)
        right = min(chromatogram.shape[0], _r + window / 4)
        return right
    @staticmethod
    def __resolve_boarder_overlap(boarders: np.ndarray, peak_indices: np.ndarray, signal:np.ndarray) -> tuple[
        ndarray, list[int]|None]:

        '''
        Returns the boarders of the peaks in a chromatogram.

        Args:
            boarders (np.ndarray): Boarders of the peaks.

        Returns:
            np.ndarray: Boarders of the peaks.
        '''
        new_boarders = copy(boarders)
        flag_peaks = []
        for i, boarder in enumerate(boarders):
            if i > 0:
                if boarder[0] < boarders[i - 1][1]:
                    window = signal[peak_indices[i-1]:peak_indices[i]]
                    smooth_window = gaussian_filter1d(window, 10)
                    minima = argrelmin(smooth_window, order=10)[0]
                    if not len(minima) > 0:
                        continue
                    else:
                        flag_peaks += [i-1, i]
                        y_values = np.take(smooth_window, minima)
                        new_boarder = minima[np.argmin(y_values)] + peak_indices[i-1]
                        new_boarders[i][0] = new_boarder
                        new_boarders[i-1][1] = new_boarder
        if flag_peaks:
            flag_peaks = list(set(flag_peaks))
        return new_boarders, flag_peaks

    @staticmethod
    def __calculate_areas(time: np.ndarray, signal: np.ndarray, boarders: np.ndarray) -> list[float]:

        '''
        Returns the areas of the peaks in a signal.

        Args:
            time (np.ndarray): Retention times of the signal in minutes.
            signal (np.ndarray): Signal to calculate the areas for.
            boarders (np.ndarray): Boarders of the peaks.

        Returns:
            list[float]: Areas of the peaks.
        '''

        areas = []
        for boarder in boarders:
            area = Peak_Detection_FID.peak_area(time, signal, boarder[0], boarder[1])
            areas.append(area)
        return areas
    @staticmethod
    def __initialize_peaks(peak_rts: np.ndarray, peak_heights: np.ndarray, peak_widths: np.ndarray,
                           peak_boarders: np.ndarray, peak_areas: list[float], flag_peaks_list: list[int]|None) -> dict[float, FID_Peak]:
        '''
        Returns a dictionary of FID peaks.

        Args:
            peak_rts (np.ndarray): Retention times of the peaks.
            peak_heights (np.ndarray): Heights of the peaks.
            peak_widths (np.ndarray): Widths of the peaks.
            peak_boarders (np.ndarray): Boarders of the peaks.
            peak_areas (list[float]): Areas of the peaks.
            flag_peaks_list (list[int]|None): List of indices of peaks to flag for boarder overlap.

        Returns:
            dict[float, FID_Peak]: Dictionary of FID peaks.
        '''

        peaks = {}
        for i, rt in enumerate(peak_rts):
            if flag_peaks_list:
                if i in flag_peaks_list:
                    flag =  "overlap"
                else:
                    flag = None
            else:
                flag = None
            rt = round(rt, 3)

            peak = FID_Peak(rt, peak_heights[i], peak_widths[i], np.array([peak_boarders[i][0], peak_boarders[i][1]]),
                            peak_areas[i], flag)
            peaks[peak.rt] = peak
        return peaks
