'''Chromatogram: the detector signal of one injection as an object instead of a (2, N) array.

The agent layer holds chromatograms behind handles, so they carry their own axis, units and
derived quantities, and refuse a malformed axis when they are built rather than when an algorithm
indexes it.
'''

import pickle

import numpy as np
import pytest

from pygecko.gc_tools import Chromatogram
from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.injection.ms_injection import MS_Injection

from .conftest import make_fid_injection, make_ms_injection


def make_chromatogram(signal, dt=0.01, kind='FID'):
    return Chromatogram(np.arange(len(signal)) * dt, np.asarray(signal, float), kind=kind)


class TestValidation:

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match='equal length'):
            Chromatogram([0.0, 0.1, 0.2], [1.0, 2.0], kind='FID')

    @pytest.mark.parametrize('time', [[0.0, 0.2, 0.1], [0.0, 0.1, 0.1]], ids=['decreasing', 'repeated'])
    def test_a_time_axis_that_is_not_strictly_increasing_raises(self, time):
        with pytest.raises(ValueError, match='strictly increasing'):
            Chromatogram(time, [1.0, 2.0, 3.0], kind='FID')

    def test_a_two_dimensional_axis_raises(self):
        with pytest.raises(ValueError):
            Chromatogram(np.zeros((2, 3)), np.zeros((2, 3)), kind='FID')


class TestDerivedQuantities:

    def test_scan_rate_is_the_sampling_interval_in_minutes(self):
        assert make_chromatogram([1.0] * 10, dt=0.005).scan_rate == pytest.approx(0.005)

    def test_run_time_spans_the_axis(self):
        assert make_chromatogram([1.0] * 11, dt=0.1).run_time == pytest.approx((0.0, 1.0))

    def test_a_complete_signal_has_no_empty_ranges(self):
        assert make_chromatogram([10.0] * 50).empty_ranges() == []

    def test_contiguous_zero_points_are_one_range(self):
        signal = [10.0] * 50
        signal[10:13] = [0.0, 0.0, 0.0]
        assert make_chromatogram(signal).empty_ranges() == [pytest.approx((0.10, 0.12))]

    def test_every_gap_is_reported(self):
        signal = [10.0] * 100
        for i in (10, 20, 21, 22, 40):
            signal[i] = 0.0
        assert len(make_chromatogram(signal).empty_ranges()) == 3


class TestSlicing:

    def test_a_slice_keeps_the_processed_signal_aligned(self):
        chromatogram = make_chromatogram([1.0, 2.0, 3.0, 4.0])
        chromatogram.processed = np.array([10.0, 20.0, 30.0, 40.0])
        tail = chromatogram[2:]
        np.testing.assert_allclose(tail.time, [0.02, 0.03])
        np.testing.assert_allclose(tail.intensity, [3.0, 4.0])
        np.testing.assert_allclose(tail.processed, [30.0, 40.0])
        assert tail.kind == 'FID'

    def test_a_slice_of_an_unprocessed_chromatogram_is_unprocessed(self):
        assert make_chromatogram([1.0, 2.0, 3.0, 4.0])[1:].processed is None


class TestOnTheInjection:

    def test_an_fid_injection_holds_a_chromatogram(self):
        injection = make_fid_injection()
        assert isinstance(injection.chromatogram, Chromatogram)
        assert injection.chromatogram.kind == 'FID'
        assert injection.chromatogram.processed is None

    def test_pick_peaks_fills_the_processed_signal(self):
        injection = make_fid_injection()
        injection.pick_peaks()
        assert len(injection.chromatogram.processed) == len(injection.chromatogram.time)

    def test_the_solvent_delay_crop_is_a_chromatogram_slice(self):
        injection = make_fid_injection(solvent_delay=2.0)
        assert isinstance(injection.chromatogram, Chromatogram)
        assert injection.chromatogram.time[0] == pytest.approx(2.0, abs=0.01)

    def test_the_injection_is_silent_about_missing_signal_at_construction(self, recwarn, capsys):
        injection = make_fid_injection()
        injection.chromatogram.empty_ranges()
        assert capsys.readouterr().out == ''
        assert len(recwarn) == 0


def slot_state(obj):
    '''Returns the slot values of obj as pickle stores them, keyed by name.'''
    names = [name for cls in type(obj).__mro__ for name in getattr(cls, '__slots__', ())]
    return {name: getattr(obj, name) for name in names if hasattr(obj, name)}


class TestOldPickles:
    '''Injections pickled while the chromatogram was a (2, N) array still load.'''

    def test_an_ndarray_chromatogram_is_wrapped_and_the_stale_processed_one_dropped(self):
        injection = make_fid_injection(solvent_delay=0.001)
        old_layout = slot_state(injection)
        old_layout['chromatogram'] = np.array([injection.chromatogram.time, injection.chromatogram.intensity])
        old_layout['processed_chromatogram'] = old_layout['chromatogram'][:, 10:]
        old_layout.pop('_baseline_settings', None)

        restored = object.__new__(FID_Injection)
        restored.__setstate__((None, old_layout))

        assert isinstance(restored.chromatogram, Chromatogram)
        assert restored.chromatogram.kind == 'FID'
        assert restored.chromatogram.processed is None
        restored.pick_peaks()
        assert len(restored.peaks) == 2

    def test_an_ms_injection_gets_a_tic(self):
        injection = make_ms_injection([])
        old_layout = slot_state(injection)
        old_layout['chromatogram'] = np.array([injection.chromatogram.time, injection.chromatogram.intensity])

        restored = object.__new__(MS_Injection)
        restored.__setstate__((None, old_layout))

        assert restored.chromatogram.kind == 'TIC'

    def test_a_chromatogram_round_trips_through_pickle(self):
        injection = make_fid_injection()
        injection.pick_peaks()
        restored = pickle.loads(pickle.dumps(injection))
        np.testing.assert_array_equal(restored.chromatogram.processed, injection.chromatogram.processed)
