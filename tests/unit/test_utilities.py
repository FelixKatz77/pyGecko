'''Tests for the Utilities namespace class.'''

import numpy as np
import pytest

from pygecko.gc_tools.utilities import Utilities


def make_chromatogram(signal, dt=0.1):
    '''Builds a (2, N) chromatogram with a regular time axis in minutes.'''
    signal = np.asarray(signal, dtype=float)
    return np.array([np.arange(len(signal)) * dt, signal])


class TestFindEmptyRanges:

    def test_all_positive_signal_has_no_empty_ranges(self):
        assert Utilities.find_empty_ranges(make_chromatogram([10.0] * 10)) == []

    def test_single_zero_marks_only_that_point(self):
        '''A single dropout must not mark the whole chromatogram as empty.'''
        signal = [10.0] * 10
        signal[5] = 0.0
        assert Utilities.find_empty_ranges(make_chromatogram(signal)) == [
            pytest.approx((0.5, 0.5))]

    def test_contiguous_zeros_give_one_range(self):
        signal = [10.0] * 10
        signal[4:7] = [0.0, 0.0, 0.0]
        assert Utilities.find_empty_ranges(make_chromatogram(signal)) == [
            pytest.approx((0.4, 0.6))]

    def test_two_gaps_give_two_ranges(self):
        signal = [10.0] * 10
        signal[1] = 0.0
        signal[7:9] = [0.0, 0.0]
        ranges = Utilities.find_empty_ranges(make_chromatogram(signal))
        assert ranges == [pytest.approx((0.1, 0.1)), pytest.approx((0.7, 0.8))]

    def test_nan_counts_as_empty(self):
        signal = [10.0] * 10
        signal[3] = np.nan
        assert Utilities.find_empty_ranges(make_chromatogram(signal)) == [
            pytest.approx((0.3, 0.3))]

    def test_threshold_is_applied(self):
        '''The documented threshold parameter must actually widen what counts as empty.'''
        signal = [10.0] * 10
        signal[4:6] = [2.0, 3.0]
        assert Utilities.find_empty_ranges(make_chromatogram(signal)) == []
        assert Utilities.find_empty_ranges(make_chromatogram(signal), threshold=5.0) == [
            pytest.approx((0.4, 0.5))]

    def test_min_duration_filters_short_gaps(self):
        signal = [10.0] * 10
        signal[2] = 0.0
        signal[5:9] = [0.0, 0.0, 0.0, 0.0]
        ranges = Utilities.find_empty_ranges(make_chromatogram(signal), min_duration=0.2)
        assert ranges == [pytest.approx((0.5, 0.8))]

    def test_empty_chromatogram_returns_no_ranges(self):
        assert Utilities.find_empty_ranges(np.array([[], []])) == []
