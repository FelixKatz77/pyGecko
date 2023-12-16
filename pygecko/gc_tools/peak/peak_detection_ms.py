import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from pygecko.gc_tools.analysis.analysis_settings import Analysis_Settings
from pygecko.gc_tools.peak.fid_peak import FID_Peak
from pygecko.gc_tools.peak.ms_peak import MS_Peak
from pygecko.gc_tools.utilities import Utilities


class Peak_Detection_MS:
    '''
    A class wrapping functions to detect peaks in MS chromatograms.
    '''

    @staticmethod
    def pick_peaks(chromatogram: np.ndarray, scans: pd.DataFrame, analysis_settings: Analysis_Settings) -> dict[float:MS_Peak]:

        '''
        Returns a dictionary of MS peaks.

        Args:
            chromatogram (np.ndarray): Chromatogram to detect peaks in.
            scans (pd.DataFrame): Mass traces of the chromatogram.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            dict[float:MS_Peak]: Dictionary of MS peaks.
        '''

        peak_indices, peak_rts, peak_heights = Peak_Detection_MS.__detect_peaks_scipy(chromatogram, analysis_settings)
        spectra = Peak_Detection_MS.__extract_mass_spectrum(scans, peak_rts, peak_indices, analysis_settings)
        peaks = Peak_Detection_MS.__initialize_peaks(peak_rts, peak_heights, spectra)
        return peaks

    @staticmethod
    def __detect_peaks_scipy(chromatogram: np.ndarray, analysis_settings: Analysis_Settings) -> tuple[
        np.ndarray, np.ndarray, np.ndarray]:

        '''
        Returns the peak indices, retention times and heights of a chromatogram.

        Args:
            chromatogram (np.ndarray): Chromatogram to detect peaks in.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: Peak indices, retention times and heights.
        '''

        index = chromatogram[0]
        chromatogram = chromatogram[1]

        min_height = analysis_settings.pop('height', np.min(chromatogram) * 50)
        prominence = analysis_settings.pop('prominence_ms', 1)
        prominence = np.median(chromatogram) * prominence
        width = analysis_settings.pop('width', 0)

        peak_indices, peak_properties = find_peaks(chromatogram,
                                                   prominence=prominence, width=width, height=min_height)
        peak_heights = peak_properties['peak_heights']
        peak_rts = index[peak_indices]
        return peak_indices, peak_rts, peak_heights

    @staticmethod
    def __extract_mass_spectrum(scans, peak_rts, peak_indices, analysis_settings: Analysis_Settings):

        '''
        Returns the mass spectra for the peaks of a chromatogram.

        Args:
            scans (pd.DataFrame): Mass traces of the chromatogram.
            peak_rts (np.ndarray): Retention times of the peaks.
            peak_indices (np.ndarray): Indices of the chromatogram at which peaks are located.
            analysis_settings (Analysis_Settings): Data_Method object containing settings for the peak detection.

        Returns:
            dict[float:dict[float:float]]: Mass spectra of the peaks.
        '''

        prominence = analysis_settings.pop('trace_prominence', 500)

        mass_spectrums = {}
        for i in peak_rts:
            mass_spectrums[i] = {}
        for mass_trace in scans:
            rts = scans.index.to_numpy() / 60000
            chromatogram = scans[mass_trace].to_numpy().transpose()
            indices, properties = find_peaks(chromatogram, prominence=prominence)
            for peak_index in peak_indices:
                for index in indices:
                    if Utilities.check_interval(index, peak_index, 5):
                        mass_spectrums[rts[peak_index]].update({mass_trace: chromatogram[peak_index]})
        return mass_spectrums

    @staticmethod
    def __initialize_peaks(peak_rts: np.ndarray, peak_heights:np.ndarray, mass_spectra:dict[float:dict[float:float]]) -> dict[float:FID_Peak]:

        '''
        Returns a dictionary of MS peaks.

        Args:
            peak_rts (np.ndarray): Retention times of the peaks.
            peak_heights (np.ndarray): Heights of the peaks.
            mass_spectra (dict[float:dict[float:float]]): Mass spectra of the peaks.

        Returns:
            dict[float:FID_Peak]: Dictionary of MS peaks.
        '''

        peaks = {}
        for i, rt in enumerate(peak_rts):
            rt_min = round(rt, 3)
            intensities = list(mass_spectra[rt].values())
            rel_intensities = np.divide(intensities, np.max(intensities)) * 100
            l = [(i, j, k) for i, j, k in zip(mass_spectra[rt].keys(), intensities, rel_intensities)]
            mass_spectrum = np.array(l,
                                     dtype=dict(names=['mz', 'intensity', 'rel_intensity'], formats=['f8', 'f8', 'f8']))
            peak = MS_Peak(rt_min, peak_heights[i], mass_spectrum)
            peaks[peak.rt] = peak
        return peaks
