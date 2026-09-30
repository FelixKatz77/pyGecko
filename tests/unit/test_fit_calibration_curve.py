'''Analysis.fit_calibration_curve returns the fit quality with the fit.'''

import pytest

from pygecko.analysis import Analysis
from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.sequence.fid_sequence import FID_Sequence

from .conftest import make_injection, make_peak

IS_RT = 4.0
ANALYTE_RT = 6.0


def make_calibration_sequence(analyte_areas):
    '''Builds a sequence with one injection per analyte area, each against a standard of area 1.'''
    injections = {}
    for i, area in enumerate(analyte_areas):
        name = f'CAL-{i}'
        injection = make_injection([make_peak(IS_RT, height=1.0), make_peak(ANALYTE_RT, height=area)],
                                   sample_name=name)
        injection.set_internal_standard(IS_RT, name='IS', smiles='CCCCCCCCCCCC')
        injections[name] = injection
    return FID_Sequence({}, injections)


def test_returns_slope_intercept_and_r_squared():
    sequence = make_calibration_sequence([1.0, 2.0, 3.0])
    slope, intercept, r2 = Analysis.fit_calibration_curve(sequence, Analyte(ANALYTE_RT), [0.5, 1.0, 1.5])
    assert slope == pytest.approx(0.5)
    assert intercept == 0.0
    assert r2 == pytest.approx(1.0)


def test_a_poor_fit_warns():
    sequence = make_calibration_sequence([1.0, 2.0, 3.0])
    with pytest.warns(UserWarning, match='not accurate'):
        Analysis.fit_calibration_curve(sequence, Analyte(ANALYTE_RT), [1.5, 0.1, 1.0])


def test_the_standard_keeps_its_analyte():
    # flag_peak(internal_standard.rt) used to overwrite the standard's analyte with None, which made
    # a following polyarc quantification raise 'No molecule assigned for Standard'.
    sequence = make_calibration_sequence([1.0, 2.0, 3.0])
    Analysis.fit_calibration_curve(sequence, Analyte(ANALYTE_RT), [0.5, 1.0, 1.5])
    for injection in sequence.injections.values():
        assert injection[IS_RT].analyte is injection.internal_standard
