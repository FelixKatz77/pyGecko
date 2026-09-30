'''Quantification returns yields as float percentages; rounding is the caller's decision.

int(round(yield_, 0)) turned 0.4 % into 0 and 99.6 % into 100 before any caller could see it.
'''

import pytest

from pygecko.gc_tools.analysis.quantification import Quantification
from pygecko.gc_tools.analyte import Analyte

from .conftest import make_peak


def peak_with(area, smiles=None, rt=5.0):
    peak = make_peak(rt, height=area)
    if smiles:
        peak.analyte = Analyte(rt, smiles=smiles)
    return peak


def test_ratio_keeps_the_fraction():
    yield_ = Quantification.quantify_ratio(peak_with(1.0), peak_with(3.0))
    assert isinstance(yield_, float)
    assert yield_ == pytest.approx(100 / 3)


def test_a_small_ratio_is_not_rounded_to_zero():
    assert Quantification.quantify_ratio(peak_with(0.004), peak_with(1.0)) == pytest.approx(0.4)


def test_calibration_keeps_the_fraction():
    yield_ = Quantification.quantify_calibration(peak_with(1.0), peak_with(3.0), slope=1.0, intercept=0.0)
    assert yield_ == pytest.approx(100 / 3)


def test_polyarc_keeps_the_fraction():
    # Area per carbon: 1/3 for the C3 analyte against 3/6 for the C6 standard.
    yield_ = Quantification.quantify_polyarc(peak_with(1.0, 'CCC'), peak_with(3.0, 'CCCCCC'))
    assert yield_ == pytest.approx(200 / 3)


@pytest.mark.parametrize('analyte_smiles, standard_smiles, missing', [
    (None, 'CCCCCC', 'Analyte'),
    ('CCC', None, 'Standard'),
])
def test_polyarc_needs_a_molecule_on_both_peaks(analyte_smiles, standard_smiles, missing):
    analyte, standard = peak_with(1.0), peak_with(3.0)
    analyte.analyte = Analyte(5.0, smiles=analyte_smiles)
    standard.analyte = Analyte(4.0, smiles=standard_smiles)
    with pytest.raises(TypeError, match=f'No molecule assigned for {missing}'):
        Quantification.quantify_polyarc(analyte, standard)
