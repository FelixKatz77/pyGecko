'''FID_Injection.quantify dispatches every method Quantification offers and fails loudly otherwise.

It used to handle polyarc and calibration only and return None for anything else, including
'ratio', so a typo or an unsupported method reached the plate as a silent missing yield.
'''

import pytest

from .conftest import make_fid_injection

IS_RT = 4.0
ANALYTE_RT = 6.0


@pytest.fixture
def picked():
    '''An FID injection with both synthetic peaks picked and the earlier one set as the standard.'''
    injection = make_fid_injection(solvent_delay=0.001)
    injection.pick_peaks()
    injection.set_internal_standard(IS_RT)
    return injection


def analyte_rt(injection):
    return min(injection.peaks, key=lambda rt: abs(rt - ANALYTE_RT))


class TestQuantifyDispatch:

    def test_ratio_returns_a_float_yield(self, picked):
        yield_ = picked.quantify(analyte_rt(picked), method='ratio')
        assert isinstance(yield_, float)
        assert yield_ == pytest.approx(100, abs=1)

    def test_calibration_applies_slope_and_intercept(self, picked):
        rt = analyte_rt(picked)
        ratio = picked.quantify(rt, method='ratio')
        assert picked.quantify(rt, method='calibration', slope=2.0, intercept=0.0) == pytest.approx(2 * ratio)

    def test_an_unknown_method_raises(self, picked):
        with pytest.raises(ValueError, match='bogus'):
            picked.quantify(analyte_rt(picked), method='bogus')

    def test_a_missing_internal_standard_raises(self):
        injection = make_fid_injection(solvent_delay=0.001)
        injection.pick_peaks()
        with pytest.raises(ValueError, match='no internal standard'):
            injection.quantify(analyte_rt(injection), method='ratio')

    def test_an_unknown_retention_time_raises_a_descriptive_key_error(self, picked):
        with pytest.raises(KeyError, match='no peak at 6.5'):
            picked.quantify(6.5, method='ratio')
