import numpy as np
from pybaselines import Baseline
from scipy.integrate import simpson
from scipy.signal import find_peaks, savgol_filter, argrelmin
from scipy.ndimage import gaussian_filter1d
from statsmodels.stats.stattools import durbin_watson
from copy import copy
from pygecko.gc_tools.analysis.analysis_settings import Analysis_Settings
from pygecko.gc_tools.peak.fid_peak import FID_Peak



class Peak_Detection_FID:
    '''
    A class wrapping functions to detect peaks in FID chromatograms.
    '''

    @staticmethod
    def baseline_correction(chromatogram: np.ndarray, analysis_settings: Analysis_Settings) -> np.ndarray:

        '''
        Returns a baseline corrected chromatogram.

        Args:
            chromatogram (np.ndarray): Chromatogram to correct.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the baseline correction.

        Returns:
            np.ndarray: Baseline corrected chromatogram.
        '''

        indices_range = analysis_settings.indices_range
        chrom_smooth = Peak_Detection_FID.__savgol(chromatogram[:, indices_range[0]:indices_range[1]],
                                                   analysis_settings)
        chrom_corr, baseline = Peak_Detection_FID.__baseline_filter(chrom_smooth, analysis_settings)
        return chrom_corr

    @staticmethod
    def pick_peaks(chromatogram: np.ndarray, analysis_settings: Analysis_Settings) -> dict[float:FID_Peak]:

        '''
        Returns a dictionary of FID peaks.

        Args:
            chromatogram (np.ndarray): Chromatogram to detect peaks in.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            dict[float:FID_Peak]: Dictionary of FID peaks.
        '''

        peak_rts, peak_widths, peak_heights, peak_boarders, peak_areas = Peak_Detection_FID.__detect_peaks(chromatogram,
                                                                                                           analysis_settings)
        peaks = Peak_Detection_FID.__initialize_peaks(peak_rts, peak_heights, peak_widths, peak_boarders, peak_areas)
        return peaks

    @staticmethod
    def __detect_peaks(chrom_corr: np.ndarray, analysis_settings: Analysis_Settings):

        '''
        Returns the peak retention times, widths, heights, boarders and areas of a chromatogram.

        Args:
            chrom_corr (np.ndarray): Chromatogram to detect peaks in.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[float]]: Peak retention times, widths, heights,
            boarders and areas.
        '''

        indices_range = analysis_settings.pop('indices_range', [None, None])
        prominence = analysis_settings.pop('prominence_fid', np.mean(chrom_corr[1]))
        width = analysis_settings.pop('width', 0)
        #TODO: Implement S/N.
        height = analysis_settings.pop('height', 0)

        peak_indices, peak_properties = find_peaks(chrom_corr[1], prominence=prominence, width=width,
                                                   height=height)
        peak_boarders = Peak_Detection_FID.__detect_boarders(chrom_corr, peak_indices, peak_properties['widths'],
                                                             analysis_settings)
        peak_boarders = Peak_Detection_FID.__resolve_boarder_overlap(peak_boarders, peak_indices, chrom_corr)
        peak_areas = Peak_Detection_FID.__calculate_areas(chrom_corr, peak_boarders)
        peak_widths = [((boarder[1] - boarder[0])*analysis_settings.scan_rate) for boarder in peak_boarders]
        peak_indices = peak_indices + indices_range[0]
        peak_boarders = [((boarder + indices_range[0])*analysis_settings.scan_rate) + chrom_corr[0][0] for boarder in peak_boarders]
        peak_rts = chrom_corr[0][peak_indices]
        return peak_rts, peak_widths, peak_properties['peak_heights'], peak_boarders, peak_areas

    @staticmethod
    def __baseline_filter(chromatogram: np.ndarray, analysis_settings: Analysis_Settings) -> tuple[np.ndarray, np.ndarray]:

        '''
        Returns a baseline corrected chromatogram and the baseline.

        Args:
            chromatogram (np.ndarray): Chromatogram to correct.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the baseline correction.

        Returns:
            tuple[np.ndarray, np.ndarray]: Baseline corrected chromatogram and baseline.
        '''

        max_half_window = analysis_settings.pop('max_half_window', 200)
        baseline_fitter = Baseline(x_data=chromatogram[0])
        baseline = baseline_fitter.snip(chromatogram[1], max_half_window=max_half_window)[0]
        y_corr = chromatogram[1] - baseline
        return np.vstack((chromatogram[0], y_corr)), baseline

    @staticmethod
    def __savgol(chromatogram: np.ndarray, analysis_settings: Analysis_Settings) -> np.ndarray:

        '''
        Returns a Savitzky-Golay filtered chromatogram.

        Args:
            chromatogram (np.ndarray): Chromatogram to filter.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the Savitzky-Golay filter.

        Returns:
            np.ndarray: Savitzky-Golay filtered chromatogram.
        '''

        savgol_window = analysis_settings.pop('savgol_window', Peak_Detection_FID.__optimize_savgol_window(chromatogram))
        y_smooth = savgol_filter(chromatogram[1], savgol_window, 2)
        return np.vstack((chromatogram[0], y_smooth))

    @staticmethod
    def __optimize_savgol_window(chromatogram: np.ndarray) -> int:

        '''
        Returns the optimal Savitzky-Golay window size for a chromatogram identified using the Durbin-Watson statistic.

        Args:
            chromatogram: Chromatogram to optimize the Savitzky-Golay window size for.

        Returns:
            int: Optimal Savitzky-Golay window size.
        '''

        windows = [5, 7, 9, 11, 21, 31]
        best_dw = 0
        for window in windows:
            _y = savgol_filter(chromatogram[1], window, 2)
            resids = chromatogram[1] - _y
            dw = durbin_watson(resids)
            if abs(2 - dw) < abs(2 - best_dw):
                best_dw = dw
                best_window = window
        return best_window

    @staticmethod
    def __detect_boarders(chromatogram: np.ndarray, peak_indices: np.ndarray, peak_widths,
                          analysis_settings: Analysis_Settings) -> np.ndarray:

        '''
        Returns the boarders of the peaks in a chromatogram.

        Args:
            chromatogram (np.ndarray): Chromatogram to detect the boarders for.
            peak_indices (np.ndarray): Indices of the chromatogram at which peaks are located.
            peak_widths (np.ndarray): Widths of the peaks.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the boarder detection.

        Returns:
            np.ndarray: Boarders of the peaks.
        '''

        first_diff = np.diff(chromatogram[1]) / analysis_settings.scan_rate
        threshold = analysis_settings.pop('boarder_threshold', abs(np.mean(first_diff))*0.5)
        window = analysis_settings.pop('boarder_window', 100)
        boarders = np.empty((peak_indices.shape[0], 2), int)
        for i, index in enumerate(peak_indices):
            _l = int(round(index - peak_widths[i] / 2, 0))
            interval = first_diff[_l - window:_l]
            while np.mean(interval) > threshold:
                _l -= 1
                interval = first_diff[_l - window:_l]
            left = _l - window / 4
            if left < 0:
                left = 0
            _r = int(round(index + peak_widths[i] / 2, 0))
            interval = first_diff[_r:_r + window]
            while abs(np.mean(interval)) > threshold:
                _r += 1
                interval = first_diff[_r:_r + window]
            right = _r + window / 4
            if right > chromatogram.shape[1]:
                right = chromatogram.shape[1]
            boarder = [left, right]
            boarders[i] = boarder
        return boarders

    @staticmethod
    def __resolve_boarder_overlap(boarders: np.ndarray, peak_indices: np.ndarray, chromatogram:np.ndarray) -> np.ndarray:

        '''
        Returns the boarders of the peaks in a chromatogram.

        Args:
            boarders (np.ndarray): Boarders of the peaks.

        Returns:
            np.ndarray: Boarders of the peaks.
        '''
        new_boarders = copy(boarders)
        for i, boarder in enumerate(boarders):
            if i > 0:
                if boarder[0] < boarders[i - 1][1]:
                    window = chromatogram[1][peak_indices[i-1]:peak_indices[i]]
                    smooth_window = gaussian_filter1d(window, 10)
                    minima = argrelmin(smooth_window, order=10)[0]
                    if not len(minima) > 0:
                        return boarders
                    y_values = np.take(smooth_window, minima)
                    new_boarder = minima[np.argmin(y_values)] + peak_indices[i-1]
                    new_boarders[i][0] = new_boarder
                    new_boarders[i-1][1] = new_boarder
        return new_boarders

    @staticmethod
    def __calculate_areas(chromatogram: np.ndarray, boarders: np.ndarray) -> list[float]:

        '''
        Returns the areas of the peaks in a chromatogram.

        Args:
            chromatogram (np.ndarray): Chromatogram to calculate the areas for.
            boarders (np.ndarray): Boarders of the peaks.

        Returns:
            list[float]: Areas of the peaks.
        '''

        areas = []
        for boarder in boarders:
            area = simpson(chromatogram[1][boarder[0]:boarder[1]])
            areas.append(area)
        return areas

    @staticmethod
    def __set_indices_range(analysis_settings: Analysis_Settings) -> list[float|None]:

        '''
        Takes a Analysis_Settings object and returns the indices range for peak detection.
        '''

        if not analysis_settings.time_range:
            indices_range = [None, None]
        else:
            indices_range = [x / analysis_settings.scan_rate for x in analysis_settings.time_range]
        return indices_range

    @staticmethod
    def __initialize_peaks(peak_rts: np.ndarray, peak_heights: np.ndarray, peak_widths: np.ndarray,
                           peak_boarders: np.ndarray, peak_areas: list[float]) -> dict[float, FID_Peak]:
        '''
        Returns a dictionary of FID peaks.

        Args:
            peak_rts (np.ndarray): Retention times of the peaks.
            peak_heights (np.ndarray): Heights of the peaks.
            peak_widths (np.ndarray): Widths of the peaks.
            peak_boarders (np.ndarray): Boarders of the peaks.
            peak_areas (list[float]): Areas of the peaks.

        Returns:
            dict[float, FID_Peak]: Dictionary of FID peaks.
        '''

        peaks = {}
        for i, rt in enumerate(peak_rts):
            rt = round(rt, 3)
            peak = FID_Peak(rt, peak_heights[i], peak_widths[i], np.array([peak_boarders[i][0], peak_boarders[i][1]]),
                            peak_areas[i])
            peaks[peak.rt] = peak
        return peaks
