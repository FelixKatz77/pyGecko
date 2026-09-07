'''Characterisation tests pinning the MS analyte-detection thresholds.

The split-GC work loosened both thresholds globally (isotope tolerance 0.055 -> 0.15 and the
parent-ion m/z floor 2/3 -> 1/3 of the spectrum's maximum m/z) to make its own data match, at the
cost of more false-positive assignments for every other user. These tests pin the strict values as
the defaults and prove the looser split-GC values remain reachable per call, so the trade-off is an
explicit opt-in rather than a silent change.
'''

import pytest

from .conftest import make_ms_injection

BENZENE = 'c1ccccc1'
BENZENE_MZ = 78.0
# Theoretical M+1/M intensity ratio for C6H6, as computed by brainpy in __isotopic_ratio_check.
BENZENE_THEO_RATIO = 0.0655844491155197


class TestParentIonMzFloor:
    '''The parent ion must sit above a fraction of the spectrum's maximum m/z to be trusted.'''

    @pytest.fixture
    def injection(self, ms_peak_factory):
        # m/z 78 is 39% of the maximum m/z (200): above a 1/3 floor, below a 2/3 floor.
        # A fragment ion is listed first so the parent is not at index 0 (see
        # test_parent_ion_at_index_zero_passes_isotope_check).
        peak = ms_peak_factory(5.0, {51.0: 100.0, BENZENE_MZ: 1000.0, 200.0: 50.0})
        return make_ms_injection([peak])

    def test_default_floor_rejects_low_parent_ion(self, injection):
        assert injection.match_mol(BENZENE, check_iso=False) is None

    def test_loosened_floor_accepts_low_parent_ion(self, injection):
        match = injection.match_mol(BENZENE, check_iso=False, min_mz_fraction=1 / 3)
        assert match is not None
        assert match.rt == 5.0


class TestIsotopeRatioTolerance:
    '''The measured M+1/M ratio must sit within max_isotopic_diff of the theoretical ratio.'''

    @pytest.fixture
    def injection(self, ms_peak_factory):
        # Measured ratio deviates from theory by 0.10: outside a 0.055 window, inside a 0.15 one.
        m_plus_one = 1000.0 * (BENZENE_THEO_RATIO + 0.10)
        peak = ms_peak_factory(5.0, {51.0: 100.0, BENZENE_MZ: 1000.0, 79.0: m_plus_one})
        return make_ms_injection([peak])

    def test_default_tolerance_rejects_deviating_isotope_ratio(self, injection):
        assert injection.match_mol(BENZENE) is None

    def test_loosened_tolerance_accepts_deviating_isotope_ratio(self, injection):
        match = injection.match_mol(BENZENE, max_isotopic_diff=0.15)
        assert match is not None
        assert match.rt == 5.0

    def test_ratio_close_to_theory_matches_under_default(self, ms_peak_factory):
        m_plus_one = 1000.0 * (BENZENE_THEO_RATIO + 0.01)
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {51.0: 100.0, BENZENE_MZ: 1000.0, 79.0: m_plus_one})])
        assert injection.match_mol(BENZENE) is not None

    def test_parent_ion_at_index_zero_passes_isotope_check(self, ms_peak_factory):
        m_plus_one = 1000.0 * (BENZENE_THEO_RATIO + 0.01)
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {BENZENE_MZ: 1000.0, 79.0: m_plus_one})])
        assert injection.match_mol(BENZENE) is not None


class TestMatchMzHonoursSettings:
    '''match_mz applies the same parent-ion filter as match_mol and must read the same settings.

    match_mz calls analysis_settings.update(**kwargs), so a caller passing min_mz_fraction or
    min_rel_intensity has them validated and stored on the injection. If the filter then ignores them
    and uses hardcoded values, the kwarg silently does nothing where it was passed while still taking
    effect on a later match_mol against the same injection, because settings persist.
    '''

    def test_default_mz_floor_rejects_low_parent_ion(self, ms_peak_factory):
        # m/z 78 is 39% of the maximum m/z (200): above a 1/3 floor, below the default 2/3 floor.
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {51.0: 100.0, BENZENE_MZ: 1000.0, 200.0: 50.0})])
        assert injection.match_mz(BENZENE_MZ) is None

    def test_min_mz_fraction_is_honoured(self, ms_peak_factory):
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {51.0: 100.0, BENZENE_MZ: 1000.0, 200.0: 50.0})])
        match = injection.match_mz(BENZENE_MZ, min_mz_fraction=1 / 3)
        assert match is not None
        assert match.rt == 5.0

    def test_default_relative_intensity_floor_rejects_weak_ion(self, ms_peak_factory):
        # m/z 78 sits at 3% relative intensity, below the default min_rel_intensity of 4.
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {51.0: 1000.0, BENZENE_MZ: 30.0})])
        assert injection.match_mz(BENZENE_MZ) is None

    def test_min_rel_intensity_is_honoured(self, ms_peak_factory):
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {51.0: 1000.0, BENZENE_MZ: 30.0})])
        match = injection.match_mz(BENZENE_MZ, min_rel_intensity=2.0)
        assert match is not None
        assert match.rt == 5.0

    def test_strong_high_mass_ion_matches_under_defaults(self, ms_peak_factory):
        injection = make_ms_injection(
            [ms_peak_factory(5.0, {51.0: 100.0, BENZENE_MZ: 1000.0})])
        assert injection.match_mz(BENZENE_MZ) is not None
