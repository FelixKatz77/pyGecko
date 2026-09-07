'''Tests for retention-time based FID/MS peak matching (split-GC workflow).'''

import pytest

from pygecko.analysis import Analysis

from .conftest import make_injection, make_peak


class TestMatchRt:

    def test_returns_nearest_peak_within_tolerance(self):
        injection = make_injection([make_peak(5.00), make_peak(5.01), make_peak(6.00)])
        match = injection.match_rt(5.011, tolerance=0.05)
        assert match.rt == 5.01

    def test_returns_none_outside_tolerance(self):
        injection = make_injection([make_peak(5.00)])
        assert injection.match_rt(5.50, tolerance=1 / 60) is None

    def test_default_tolerance_is_one_second(self):
        injection = make_injection([make_peak(5.00)])
        assert injection.match_rt(5.0 + 0.5 / 60) is not None
        assert injection.match_rt(5.0 + 2.0 / 60) is None

    def test_applies_rt_mapping_function(self):
        # FID peak elutes 0.1 min after the MS peak; without the mapping it is out of tolerance.
        injection = make_injection([make_peak(5.10)])
        assert injection.match_rt(5.00, tolerance=1 / 60) is None
        match = injection.match_rt(5.00, func=Analysis.constant_offset(0.1), tolerance=1 / 60)
        assert match.rt == 5.10

    def test_excludes_internal_standard_by_default(self):
        injection = make_injection([make_peak(5.00, flags=['standard'])])
        assert injection.match_rt(5.00, tolerance=0.05) is None
        assert injection.match_rt(5.00, tolerance=0.05, exclude_standard=False).rt == 5.00

    def test_return_candidates_gives_list_of_peaks(self):
        injection = make_injection([make_peak(5.00), make_peak(5.02)])
        candidates = injection.match_rt(5.005, tolerance=0.05, return_candidates=True)
        assert isinstance(candidates, list)
        assert sorted(peak.rt for peak in candidates) == [5.00, 5.02]

    def test_equidistant_candidates_are_both_returned(self):
        injection = make_injection([make_peak(5.00), make_peak(5.01)])
        candidates = injection.match_rt(5.005, tolerance=0.05, return_candidates=True)
        assert sorted(peak.rt for peak in candidates) == [5.00, 5.01]

    def test_assigns_analyte_to_match(self):
        injection = make_injection([make_peak(5.00)])
        match = injection.match_rt(5.00, tolerance=0.05, analyte='c1ccccc1')
        assert match.analyte == 'c1ccccc1'


class TestRtMappingFunctions:

    @pytest.mark.parametrize('offset,rt,expected', [
        (0.0, 5.0, 5.0),
        (0.1, 5.0, 5.1),
        (-0.1, 5.0, 4.9),
    ])
    def test_constant_offset(self, offset, rt, expected):
        assert Analysis.constant_offset(offset)(rt) == pytest.approx(expected)

    def test_constant_offset_defaults_to_identity(self):
        assert Analysis.constant_offset()(7.25) == pytest.approx(7.25)

    @pytest.mark.parametrize('a,b,rt,expected', [
        (1.0, 0.0, 5.0, 5.0),
        (1.01, 0.05, 5.0, 5.1),
        (0.99, -0.05, 10.0, 9.85),
    ])
    def test_linear_drift(self, a, b, rt, expected):
        assert Analysis.linear_drift(a, b)(rt) == pytest.approx(expected)
