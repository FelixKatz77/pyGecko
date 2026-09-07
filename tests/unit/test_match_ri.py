'''Tests for retention-index based FID/MS peak matching (two-instrument workflow).'''

from pygecko.analysis import Analysis
from pygecko.gc_tools.analyte import Analyte

from .conftest import make_injection, make_peak


def make_ri_peak(rt, ri, height=100.0, flags=None):
    '''Builds a Peak carrying a retention index, which make_peak leaves as None.'''
    peak = make_peak(rt, height=height, flags=flags)
    peak.ri = ri
    return peak


class TestMatchRi:

    def test_returns_nearest_peak_within_tolerance(self):
        injection = make_injection([make_ri_peak(5.00, 1000), make_ri_peak(5.10, 1015),
                                    make_ri_peak(6.00, 1200)])
        assert injection.match_ri(1012, tolerance=20).ri == 1015

    def test_returns_none_outside_tolerance(self):
        injection = make_injection([make_ri_peak(5.00, 1000)])
        assert injection.match_ri(1100, tolerance=20) is None

    def test_ignores_peaks_without_a_retention_index(self):
        injection = make_injection([make_peak(5.00)])
        assert injection.match_ri(1000, tolerance=20) is None

    def test_return_candidates_gives_list_of_peaks(self):
        injection = make_injection([make_ri_peak(5.00, 1000), make_ri_peak(5.10, 1015)])
        candidates = injection.match_ri(1005, tolerance=20, return_candidates=True)
        assert isinstance(candidates, list)
        assert sorted(peak.ri for peak in candidates) == [1000, 1015]

    def test_equidistant_candidates_are_both_returned(self):
        injection = make_injection([make_ri_peak(5.00, 1000), make_ri_peak(5.10, 1010)])
        candidates = injection.match_ri(1005, tolerance=20, return_candidates=True)
        assert sorted(peak.ri for peak in candidates) == [1000, 1010]

    def test_return_candidates_keeps_the_internal_standard(self):
        '''Analysis relies on the candidate list still containing the standard; only the
        single-match branch drops it.'''
        injection = make_injection([make_ri_peak(5.00, 1000, flags=['standard']),
                                    make_ri_peak(5.10, 1010)])
        candidates = injection.match_ri(1005, tolerance=20, return_candidates=True)
        assert sorted(peak.ri for peak in candidates) == [1000, 1010]

    def test_single_match_drops_the_internal_standard(self):
        injection = make_injection([make_ri_peak(5.00, 1000, flags=['standard']),
                                    make_ri_peak(5.10, 1010)])
        assert injection.match_ri(1002, tolerance=20).ri == 1010

    def test_assigns_analyte_to_match(self):
        injection = make_injection([make_ri_peak(5.00, 1000)])
        assert injection.match_ri(1000, tolerance=20, analyte='c1ccccc1').analyte == 'c1ccccc1'


class TestFindBestRiMatch:
    '''Analysis.__find_best_ri_match disambiguates a candidate list by height ratio to the
    internal standard, so it must see every candidate match_ri/match_rt found.'''

    @staticmethod
    def _find_best(candidates, injection, ms_height_ratio):
        return Analysis._Analysis__find_best_ri_match(candidates, injection, ms_height_ratio)

    def test_picks_the_candidate_matching_the_ms_height_ratio(self):
        standard = make_peak(4.00, height=100.0, flags=['standard'])
        low = make_peak(5.00, height=50.0)
        high = make_peak(5.01, height=200.0)
        injection = make_injection([standard, low, high])
        injection.internal_standard = Analyte(standard.rt, name='Dodecane')

        # An MS height ratio of 2.0 should select the peak whose FID ratio is 200/100 = 2.0.
        assert self._find_best([low, high], injection, 2.0) is high
        # ...and a ratio of 0.5 the peak at 50/100 = 0.5.
        assert self._find_best([low, high], injection, 0.5) is low

    def test_weighs_equidistant_candidates(self):
        '''The pair that previously collided on one deviation key must both be weighable.'''
        standard = make_peak(4.00, height=100.0, flags=['standard'])
        below = make_peak(5.00, height=50.0)
        above = make_peak(5.01, height=200.0)
        injection = make_injection([standard, below, above])
        injection.internal_standard = Analyte(standard.rt, name='Dodecane')

        candidates = injection.match_rt(5.005, tolerance=0.05, return_candidates=True,
                                        exclude_standard=True)
        assert len(candidates) == 2
        assert self._find_best(candidates, injection, 2.0) is above

    def test_assigns_analyte_to_the_winner(self):
        standard = make_peak(4.00, height=100.0, flags=['standard'])
        peak = make_peak(5.00, height=100.0)
        injection = make_injection([standard, peak])
        injection.internal_standard = Analyte(standard.rt, name='Dodecane')

        best = Analysis._Analysis__find_best_ri_match([peak], injection, 1.0, analyte='c1ccccc1')
        assert best.analyte == 'c1ccccc1'
