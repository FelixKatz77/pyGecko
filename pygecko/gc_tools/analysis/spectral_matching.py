import numpy as np
from scipy.spatial import distance
from typing import TYPE_CHECKING, TypeVar


if TYPE_CHECKING:
    from pygecko.gc_tools.peak import MS_Peak
T = TypeVar('T', bound='Specral_Match')


class Spectral_Match:

    '''
    A class to represent a match between two mass spectra.

    Attributes:
        peaks (tuple[MS_Peak, MS_Peak]): Tuple containing the matched peaks.
        ms_score (float): Score of the match based on the cosine similarity of the mass spectra.
        rt_score (float): Score of the match based on the difference in retention time.
    '''

    peaks: tuple['MS_Peak', 'MS_Peak']
    ms_score: float
    rt_score: float

    __slots__ = 'peaks', 'ms_score', 'rt_score'

    def __init__(self, peaks:tuple['MS_Peak', 'MS_Peak'], ms_score, rt_score):
        self.peaks = peaks
        self.ms_score = ms_score
        self.rt_score = rt_score

    @classmethod
    def match_peaks(cls, peak1: 'MS_Peak', peak2: 'MS_Peak', rt_threshold=0.1, ms_threshold=0.9) -> T|None:

        '''
        Returns a Spectral_Match object if the peaks match and None otherwise.

        Args:
            peak1 (MS_Peak): First peak to compare.
            peak2 (MS_Peak): Second peak to compare.
            rt_threshold (float, optional): Retention time threshold for the match. Defaults to 0.1.
            ms_threshold (float, optional): Mass spectrum similarity threshold for the match. Defaults to 0.9.

        Returns:
            Spectral_Match|None: Spectral_Match object if the peaks match and None otherwise.
        '''

        rt_score = cls.__get_rt_score((peak1.rt, peak2.rt), rt_threshold)
        if rt_score > 0:
            ms_score = cls.__get_ms_score((peak1.mass_spectrum, peak2.mass_spectrum))
            if ms_score >= ms_threshold:
                return cls((peak1, peak2), ms_score, rt_score)
            else:
                return None
        else:
            return None

    @staticmethod
    def __get_rt_score(peak_rts: tuple[float, float], rt_threshold: float) -> float:

        '''
        Returns the retention time score for the comparison of two retention times.

        Args:
            peak_rts (tuple[float, float]): Tuple containing the retention times.
            rt_threshold (float): Retention time threshold for the match.

        Returns:
            float: Retention time score for the comparison of the two retention times.
        '''

        rt_diff = abs(peak_rts[0] - peak_rts[1])
        rt_score = 1 - rt_diff / rt_threshold
        return rt_score

    @staticmethod
    def __get_ms_score(ms_specs: tuple[np.ndarray, np.ndarray]) -> float:

        '''
        Returns the mass spectrum score for the comparison of two mass spectra.

        Args:
            ms_specs (tuple[np.ndarray, np.ndarray]): Tuple containing the mass spectra.

        Returns:
            float: Mass spectrum score for the comparison of the two mass spectra.
        '''

        ms_score = Spectral_Match.__get_cosine_similarity(ms_specs[0], ms_specs[1])
        return ms_score

    @staticmethod
    def __get_cosine_similarity(ms_spectrum1: np.ndarray, ms_spectrum2: np.ndarray) -> float:

        '''
        Returns the cosine similarity of two mass spectra.

        Args:
            ms_spectrum1 (np.ndarray): First mass spectrum.
            ms_spectrum2 (np.ndarray): Second mass spectrum.

        Returns:
            float: Cosine similarity of the two mass spectra.
        '''

        peak_weighted_intensities_vector, control_weighted_intensities_vector = Spectral_Match.__get_weighted_intensities_vectors(
            ms_spectrum1, ms_spectrum2)
        weighted_cosine_correlation = Spectral_Match.__calculate_weighted_cosine_correlation(
            peak_weighted_intensities_vector, control_weighted_intensities_vector)
        cosine_similarity = 1 - weighted_cosine_correlation
        return cosine_similarity

    @staticmethod
    def __get_weighted_intensities_vectors(ms_spectrum1: np.ndarray, ms_spectrum2: np.ndarray) -> tuple[
        np.ndarray, np.ndarray]:

        '''
        Returns the weighted intensity vectors for two mass spectra.

        Args:
            ms_spectrum1 (np.ndarray): First mass spectrum.
            ms_spectrum2 (np.ndarray): Second mass spectrum.

        Returns:
            tuple[np.ndarray, np.ndarray]: Tuple containing the weighted intensity vectors for the two mass spectra.
        '''

        max_mz = max(max(ms_spectrum1['mz']), max(ms_spectrum2['mz']))
        weighted_intensity_vector1 = np.zeros(int(max_mz + 1))
        weighted_intensity_vector2 = np.zeros(int(max_mz + 1))
        for mz in range(int(max_mz + 1)):
            weighted_intensity_vector1[mz] = Spectral_Match.__get_weighted_intensity(mz, ms_spectrum1)
            weighted_intensity_vector2[mz] = Spectral_Match.__get_weighted_intensity(mz, ms_spectrum2)
        return weighted_intensity_vector1, weighted_intensity_vector2

    @staticmethod
    def __get_weighted_intensity(mz: int, ms_spectrum: np.ndarray) -> float:

        '''
        Returns the weighted intensity of a mass spectrum at a given m/z.

        Args:
            mz (int): m/z to get the weighted intensity for.
            ms_spectrum (np.ndarray): Mass spectrum.

        Returns:
            float: Weighted intensity of the mass spectrum at the given m/z.
        '''

        if mz in ms_spectrum['mz']:
            rel_intensity = ms_spectrum[ms_spectrum['mz'] == mz]['rel_intensity'][0]
            weighted_intensity = Spectral_Match.__calculate_weighted_intensity(mz, rel_intensity)
        else:
            weighted_intensity = 0.0
        return weighted_intensity

    @staticmethod
    def __calculate_weighted_intensity(mz: int, rel_intensity: float) -> float:

        '''
        Returns the weighted intensity given a m/z value and relative intensity.

        Args:
            mz (int): m/z value.
            rel_intensity (float): Relative intensity.

        Returns:
            float: Weighted intensity.
        '''

        weighted_intensity = (mz ** 1.1) * (rel_intensity ** 0.5)
        return weighted_intensity

    @staticmethod
    def __calculate_weighted_cosine_correlation(weighted_intensities_vector1: np.array,
                                                weighted_intensities_vector2: np.array) -> float:

        '''
        Returns the weighted cosine correlation of two weighted intensity vectors.

        Args:
            weighted_intensities_vector1 (np.array): First weighted intensity vector.
            weighted_intensities_vector2 (np.array): Second weighted intensity vector.

        Returns:
            float: Weighted cosine correlation of the two weighted intensity vectors.
        '''

        weighted_cosine_correlation = distance.cosine(weighted_intensities_vector1,
                                                      weighted_intensities_vector2)
        return weighted_cosine_correlation