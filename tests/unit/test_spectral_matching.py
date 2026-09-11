'''
Spectral_Match scores two peaks by retention time and by the weighted cosine similarity of their
mass spectra: each m/z contributes (m/z ** 1.1) * (rel_intensity ** 0.5), absent m/z contribute 0.
'''

import numpy as np
import pytest
from scipy.spatial import distance

from pygecko.gc_tools.analysis.spectral_matching import Spectral_Match
from pygecko.gc_tools.peak import MS_Peak


def spectrum(mzs, rel_intensities):
    records = [(mz, rel * 1000, rel) for mz, rel in zip(mzs, rel_intensities)]
    return np.array(records, dtype=dict(names=['mz', 'intensity', 'rel_intensity'], formats=['f8', 'f8', 'f8']))


def peak(rt, mzs, rel_intensities):
    return MS_Peak(rt, 1000.0, 0.05, np.array([rt - 0.05, rt + 0.05]), spectrum(mzs, rel_intensities))


@pytest.fixture
def peaks():
    return (peak(5.00, [41, 57, 91, 120], [20, 100, 60, 35]),
            peak(5.03, [41, 57, 91, 105], [25, 100, 55, 10]))


def weighted_cosine(spectrum1, spectrum2):
    size = int(max(spectrum1['mz'].max(), spectrum2['mz'].max())) + 1
    vectors = []
    for spec in (spectrum1, spectrum2):
        vector = np.zeros(size)
        for mz, rel in zip(spec['mz'], spec['rel_intensity']):
            vector[int(mz)] = (mz ** 1.1) * (rel ** 0.5)
        vectors.append(vector)
    return 1 - distance.cosine(*vectors)


def test_ms_score_is_the_weighted_cosine_similarity(peaks):
    match = Spectral_Match.match_peaks(*peaks, ms_threshold=0)

    expected = weighted_cosine(peaks[0].mass_spectrum, peaks[1].mass_spectrum)
    assert match.ms_score == pytest.approx(expected, abs=1e-12)
    assert 0 < match.ms_score < 1


def test_rt_score_falls_linearly_with_the_retention_time_difference(peaks):
    match = Spectral_Match.match_peaks(*peaks, rt_threshold=0.1, ms_threshold=0)

    assert match.rt_score == pytest.approx(1 - 0.03 / 0.1)
    assert tuple(match) == peaks


def test_identical_spectra_score_one():
    first = peak(5.0, [41, 57], [50, 100])
    second = peak(5.0, [41, 57], [50, 100])

    assert Spectral_Match.match_peaks(first, second).ms_score == pytest.approx(1.0)


def test_no_match_outside_the_retention_time_threshold(peaks):
    assert Spectral_Match.match_peaks(*peaks, rt_threshold=0.02, ms_threshold=0) is None


def test_no_match_below_the_ms_threshold(peaks):
    assert Spectral_Match.match_peaks(*peaks, ms_threshold=1.0) is None
