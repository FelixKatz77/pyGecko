'''time_range is a detection window, not a crop.

baseline_correction used to crop the chromatogram to time_range before smoothing and SNIP, and
FID_Injection kept that cropped signal as its processed chromatogram. Since pick_peaks only ran the
baseline correction when no processed signal existed yet, a second pick with a wider window ran on
the first crop and could not find what the crop had cut away. MS accepted time_range and ignored it.

Now the processed signal always covers the whole stored run (solvent delay to end), peaks are found
on it, and a window keeps those whose apex lies inside it.
'''

import numpy as np
import pytest

from pygecko.gc_tools import Chromatogram
from pygecko.gc_tools.injection.ms_injection import MS_Injection

from .conftest import make_fid_chromatogram, make_fid_injection, make_scans

IS_RT = 4.0
ANALYTE_RT = 6.0


def found(injection):
    return sorted(round(rt) for rt in injection.peaks)


class TestFidWindow:

    def test_a_wider_second_window_finds_what_the_first_excluded(self):
        injection = make_fid_injection(solvent_delay=0.001)
        injection.pick_peaks(time_range=(3.0, 5.0))
        assert found(injection) == [4]
        injection.pick_peaks(time_range=(3.0, 8.0))
        assert found(injection) == [4, 6]

    def test_the_processed_signal_covers_the_whole_run_after_a_ranged_pick(self):
        injection = make_fid_injection(solvent_delay=2.0)
        injection.pick_peaks(time_range=(3.0, 5.0))
        chromatogram = injection.chromatogram
        assert chromatogram.time[0] == pytest.approx(2.0, abs=0.01)
        assert len(chromatogram.processed) == len(chromatogram.time)

    def test_a_windowed_peak_is_identical_to_the_same_peak_picked_over_the_full_run(self):
        full = make_fid_injection(solvent_delay=0.001)
        full.pick_peaks()
        windowed = make_fid_injection(solvent_delay=0.001)
        windowed.pick_peaks(time_range=(5.0, 8.0))
        (rt, peak), = windowed.peaks.items()
        assert peak.area == full.peaks[rt].area
        np.testing.assert_array_equal(peak.boarders, full.peaks[rt].boarders)

    def test_changing_only_the_window_reuses_the_processed_signal(self):
        injection = make_fid_injection(solvent_delay=0.001)
        injection.pick_peaks()
        processed = injection.chromatogram.processed
        injection.pick_peaks(time_range=(3.0, 8.0))
        assert injection.chromatogram.processed is processed

    @pytest.mark.parametrize('setting', [{'savgol_window': 31}, {'max_half_window': 50}])
    def test_changing_a_baseline_setting_recomputes_the_processed_signal(self, setting):
        injection = make_fid_injection(solvent_delay=0.001)
        injection.pick_peaks()
        processed = injection.chromatogram.processed
        injection.pick_peaks(**setting)
        assert injection.chromatogram.processed is not processed


def make_ms_injection_with_tic(peak_rts=(IS_RT, ANALYTE_RT)):
    '''Builds an MS injection whose single mass trace is the synthetic FID-shaped signal.'''
    fid = make_fid_chromatogram(peak_rts=peak_rts, points=2000)
    scans = make_scans(fid.time * 60000, [{100: value} for value in fid.intensity])
    chromatogram = Chromatogram(scans.index / 60000, scans.sum(axis=1), kind='TIC')
    return MS_Injection({'SampleName': 'A1'}, chromatogram, None, scans)


class TestMsWindow:

    def test_without_a_window_both_peaks_are_found(self):
        injection = make_ms_injection_with_tic()
        injection.pick_peaks()
        assert found(injection) == [4, 6]

    def test_a_window_excludes_the_peak_outside_it(self):
        injection = make_ms_injection_with_tic()
        injection.pick_peaks(time_range=(3.0, 5.0))
        assert found(injection) == [4]

    def test_the_kept_peak_carries_its_own_spectrum(self):
        injection = make_ms_injection_with_tic()
        injection.pick_peaks(time_range=(5.0, 8.0))
        (peak,) = injection.peaks.values()
        assert list(peak.mass_spectrum['mz']) == [100]
