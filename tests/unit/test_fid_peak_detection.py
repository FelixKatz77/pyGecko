'''Peak coordinates, windowing and re-integration in the FID detection path.

Three defects sat within a few lines of each other here, and all of them are only observable through
the peak coordinates, so they are tested together.

Peak boarders leave Peak_Detection_FID in minutes (the contract Injection._check_for_peak and the
chromatogram plot rely on), but FID_Injection.integrate used them as scan indices, so it integrated a
near-empty slice and raised on real data. The scan-index-to-minute conversion added the time-range
offset twice - once explicitly, once via the already-sliced chromatogram's own start time - which is
invisible at the default full range and wrong as soon as a time_range is set. And time_range itself
was converted to indices from absolute zero while being applied to a chromatogram already truncated
at the solvent delay, so a window was silently shifted by that delay and could select the wrong peak.

The conversion is now done by looking values up on the chromatogram's own time axis rather than
reconstructing them from scan_rate. That is exact: reconstructing an index as index * scan_rate + t0
differs from the stored axis value by around 1e-15, which is enough for np.searchsorted to flip to
the neighbouring scan, and on the CSV fixture it did so for a quarter of all indices.
'''

import numpy as np
import pytest

from .conftest import make_fid_injection

IS_RT = 4.0
ANALYTE_RT = 6.0
# Smaller than half a scan, so convert_time_to_scan rounds it to 0 and the chromatogram is not
# truncated: absolute and axis-relative windows coincide, which is what the tests below that do not
# concern windowing want.
NO_TRUNCATION = 0.001
# Several hundred scans of the 0-10 min chromatogram, so a window interpreted from absolute zero
# lands somewhere entirely different from one interpreted from the chromatogram's own start.
TRUNCATING = 2.0


def picked(**kwargs):
    '''Returns an FID injection with peaks picked over the untruncated chromatogram.'''
    injection = make_fid_injection(solvent_delay=NO_TRUNCATION)
    injection.pick_peaks(**kwargs)
    return injection


def picked_truncated(**kwargs):
    '''Returns an FID injection truncated at TRUNCATING, with peaks picked.'''
    injection = make_fid_injection(solvent_delay=TRUNCATING)
    injection.pick_peaks(**kwargs)
    return injection


class TestIntegrateRecomputesTheAreasPickPeaksSet:

    def test_integrate_does_not_raise_on_a_picked_injection(self):
        picked().integrate()

    def test_integrate_reproduces_the_areas_pick_peaks_computed(self):
        # pick_peaks integrates the baseline-corrected signal between the boarders it found, so
        # re-integrating the same peaks over the same signal has to agree. Quantification divides
        # one peak area by another and a baseline offset does not cancel between peaks of different
        # width, so an integrate that used a different signal would silently shift every yield.
        injection = picked()
        expected = {rt: peak.area for rt, peak in injection.peaks.items()}
        injection.integrate()
        for rt, peak in injection.peaks.items():
            assert peak.area == pytest.approx(expected[rt], rel=1e-9)

    def test_integrate_is_idempotent(self):
        injection = picked()
        injection.integrate()
        once = {rt: peak.area for rt, peak in injection.peaks.items()}
        injection.integrate()
        for rt, peak in injection.peaks.items():
            assert peak.area == pytest.approx(once[rt], rel=1e-9)

    def test_boarders_are_values_from_the_chromatograms_time_axis(self):
        # This is what makes integrate exact. Reconstructing a boarder as index * scan_rate + t0
        # lands a scan off wherever the reconstructed float falls on the other side of the stored
        # axis value, which shifts a narrow peak's area by percent (2.8% on the CSV fixture).
        injection = make_fid_injection(solvent_delay=NO_TRUNCATION)
        injection.pick_peaks()
        time = injection.processed_chromatogram[0]
        for peak in injection.peaks.values():
            assert peak.boarders[0] in time and peak.boarders[1] in time

    def test_integrate_produces_a_positive_area_for_every_peak(self):
        injection = picked()
        injection.integrate()
        assert injection.peaks
        for peak in injection.peaks.values():
            assert peak.area > 0

    def test_integrate_reports_an_empty_peak_list(self, capsys):
        injection = make_fid_injection(solvent_delay=NO_TRUNCATION)
        injection.integrate()
        assert 'Peaks list is empty.' in capsys.readouterr().out


class TestPeakCoordinates:

    def test_both_peaks_are_found_over_the_full_chromatogram(self):
        assert len(picked().peaks) == 2

    def test_retention_times_locate_the_peaks_in_the_signal(self):
        found = sorted(picked().peaks)
        assert found[0] == pytest.approx(IS_RT, abs=0.01)
        assert found[1] == pytest.approx(ANALYTE_RT, abs=0.01)

    def test_boarders_bracket_the_retention_time(self):
        for rt, peak in picked().peaks.items():
            assert peak.boarders[0] < rt < peak.boarders[1], f'{peak.boarders} does not bracket {rt}'

    def test_boarders_stay_inside_the_chromatogram(self):
        injection = picked()
        start, end = injection.chromatogram[0][0], injection.chromatogram[0][-1]
        for peak in injection.peaks.values():
            assert start <= peak.boarders[0] < peak.boarders[1] <= end


class TestPeakCoordinatesWithinATimeRange:

    @pytest.fixture
    def windowed(self):
        '''An injection picked over a window holding only the later of the two peaks.'''
        return picked(time_range=(5.0, 8.0))

    def test_a_time_range_selects_only_the_peaks_inside_it(self, windowed):
        assert len(windowed.peaks) == 1

    def test_the_retention_time_locates_the_peak_in_the_signal(self, windowed):
        # Before the fix the offset was added to indices that already indexed the sliced array, so
        # this either raised IndexError or reported a retention time from the wrong scan.
        assert next(iter(windowed.peaks)) == pytest.approx(ANALYTE_RT, abs=0.01)

    def test_the_boarders_bracket_the_retention_time(self, windowed):
        rt, peak = next(iter(windowed.peaks.items()))
        assert peak.boarders[0] < rt < peak.boarders[1], f'{peak.boarders} does not bracket {rt}'

    def test_the_boarders_stay_inside_the_requested_window(self, windowed):
        # The double-counted offset pushed these past the end of the run entirely.
        peak = next(iter(windowed.peaks.values()))
        assert 5.0 <= peak.boarders[0] < peak.boarders[1] <= 8.0

    def test_a_time_range_does_not_move_the_peak_it_shares_with_the_full_range(self, windowed):
        full = picked()
        assert next(iter(windowed.peaks)) == pytest.approx(sorted(full.peaks)[1], abs=0.01)

    def test_integrate_still_agrees_with_pick_peaks_within_a_time_range(self, windowed):
        expected = {rt: peak.area for rt, peak in windowed.peaks.items()}
        windowed.integrate()
        for rt, peak in windowed.peaks.items():
            assert peak.area == pytest.approx(expected[rt], rel=1e-9)


class TestTimeRangeIsRelativeToTheChromatogramStart:

    def test_a_window_around_the_earlier_peak_selects_that_peak(self):
        # The strongest discriminator: converting the window from absolute zero shifts it to
        # 5.0-7.0 on a chromatogram truncated at 2.0, so it returned the *other* peak.
        injection = picked_truncated(time_range=(3.0, 5.0))
        assert len(injection.peaks) == 1
        assert next(iter(injection.peaks)) == pytest.approx(IS_RT, abs=0.01)

    def test_a_window_around_the_later_peak_selects_that_peak(self):
        injection = picked_truncated(time_range=(5.0, 8.0))
        assert len(injection.peaks) == 1
        assert next(iter(injection.peaks)) == pytest.approx(ANALYTE_RT, abs=0.01)

    def test_the_processed_chromatogram_spans_the_requested_window(self):
        injection = picked_truncated(time_range=(3.0, 5.0))
        window = injection.processed_chromatogram[0]
        assert 3.0 <= window[0] and window[-1] <= 5.0

    def test_the_boarders_stay_inside_the_requested_window(self):
        peak = next(iter(picked_truncated(time_range=(3.0, 5.0)).peaks.values()))
        assert 3.0 <= peak.boarders[0] < peak.boarders[1] <= 5.0

    def test_a_windowed_retention_time_matches_the_full_range_one(self):
        windowed = next(iter(picked_truncated(time_range=(5.0, 8.0)).peaks))
        full = sorted(picked_truncated().peaks)[1]
        assert windowed == pytest.approx(full, abs=0.01)

    def test_a_window_starting_before_the_chromatogram_is_clamped_to_its_start(self):
        # The old arithmetic could produce a negative index here, which numpy silently reinterprets
        # as slicing from the end - a wrong window with no error.
        injection = picked_truncated(time_range=(0.5, 5.0))
        assert injection.processed_chromatogram[0][0] == pytest.approx(TRUNCATING, abs=0.01)
        assert len(injection.peaks) == 1
        assert next(iter(injection.peaks)) == pytest.approx(IS_RT, abs=0.01)

    def test_no_time_range_analyses_the_whole_chromatogram(self):
        # The default path the golden retention indices in test_fid_ri_calibration.py ride on.
        injection = picked_truncated()
        assert len(injection.peaks) == 2
        assert injection.processed_chromatogram.shape[1] == injection.chromatogram.shape[1]


class TestSolventDelayHandling:

    def test_a_zero_solvent_delay_is_honoured_rather_than_auto_detected(self):
        injection = make_fid_injection(solvent_delay=0)
        assert injection.solvent_delay == 0
        assert injection.chromatogram.shape[1] == 4000

    def test_an_omitted_solvent_delay_is_auto_detected(self):
        # This branch computed a value and then raised TypeError on None / scan_rate, because the
        # truncation used the parameter rather than the resolved self.solvent_delay.
        injection = make_fid_injection(solvent_delay=None)
        assert injection.solvent_delay > 0
        assert injection.chromatogram.shape[1] < 4000
