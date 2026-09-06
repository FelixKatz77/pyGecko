'''Tests for substrate-referenced conversion and remaining-starting-material quantification.'''

import numpy as np
import pytest
import xarray as xr

from pygecko.analysis import Analysis
from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.peak import FID_Peak
from pygecko.gc_tools.sequence.gc_sequence import GC_Sequence

from .conftest import make_ms_injection

SUBSTRATE = 'c1ccccc1'
SUBSTRATE_MZ = 78.0
IS_RT = 4.0
ANALYTE_RT = 6.0
CHROMATOGRAM = np.array([[0.0, 0.1, 0.2, 0.3], [1.0, 1.0, 1.0, 1.0]])


class StubLayout:
    '''Minimal stand-in for Reaction_Array: a one-well plate returning a fixed substrate.'''

    def __init__(self, substrate=SUBSTRATE):
        self.substrate = substrate
        self.design = xr.DataArray(np.array([[substrate]]), dims=['x', 'y'],
                                   coords={'x': ['A'], 'y': [1]})

    def get_substrate(self, pos, index=0):
        return self.substrate


@pytest.fixture
def sequences(ms_peak_factory):
    '''One-well MS and FID sequences whose analyte peak is matchable by retention time.'''
    # Fragment ion listed first so the parent is not at index 0 (see test_ms_detection_defaults).
    analyte_peak = ms_peak_factory(ANALYTE_RT, {51.0: 100.0, SUBSTRATE_MZ: 1000.0, 79.0: 65.6})
    is_peak = ms_peak_factory(IS_RT, {51.0: 100.0, 120.0: 1000.0})
    ms_injection = make_ms_injection([analyte_peak, is_peak], sample_name='SMP-A1')
    ms_injection.plate_pos = 'A1'
    ms_injection.internal_standard = Analyte(IS_RT, name='IS')

    fid_injection = FID_Injection({'SampleName': 'SMP-A1'}, CHROMATOGRAM, solvent_delay=0.1)
    fid_injection.plate_pos = 'A1'
    fid_injection.peaks = {
        ANALYTE_RT: FID_Peak(ANALYTE_RT, 100.0, 0.1, np.array([5.9, 6.1]), 100.0, None),
        IS_RT: FID_Peak(IS_RT, 100.0, 0.1, np.array([3.9, 4.1]), 100.0, None),
    }
    fid_injection.internal_standard = Analyte(IS_RT, name='IS')

    return (GC_Sequence({}, {'SMP-A1': ms_injection}),
            GC_Sequence({}, {'SMP-A1': fid_injection}))


@pytest.fixture
def fixed_remaining(monkeypatch):
    '''Pins FID_Injection.quantify so the tests exercise the conversion arithmetic, not Polyarc.'''

    def _set(remaining_pct):
        monkeypatch.setattr(FID_Injection, 'quantify',
                            lambda self, rt, **kwargs: remaining_pct)

    return _set


def conversion(ms_sequence, fid_sequence, equivalents=1.0):
    result = Analysis.calc_plate_conv(ms_sequence, fid_sequence, StubLayout(),
                                      matching='rt', rt_tolerance=0.05,
                                      equivalents=equivalents)
    return result['quantity'][0][0]


def remaining(ms_sequence, fid_sequence):
    result = Analysis.calc_plate_rsm(ms_sequence, fid_sequence, StubLayout(),
                                     matching='rt', rt_tolerance=0.05)
    return result['quantity'][0][0]


class TestConversion:

    def test_full_consumption_reads_as_complete_conversion(self, sequences, fixed_remaining):
        fixed_remaining(0.0)
        assert conversion(*sequences) == pytest.approx(100.0)

    def test_untouched_substrate_reads_as_zero_conversion(self, sequences, fixed_remaining):
        fixed_remaining(100.0)
        assert conversion(*sequences) == pytest.approx(0.0)

    def test_excess_substrate_is_floored_at_zero_not_negative(self, sequences, fixed_remaining):
        # 1.5 equiv substrate quantified against a 1 equiv standard reads ~150% remaining.
        fixed_remaining(150.0)
        assert conversion(*sequences) == pytest.approx(0.0)

    def test_equivalents_references_conversion_to_actual_loading(self, sequences, fixed_remaining):
        # 75% remaining of a 1.5 equiv charge is half the substrate consumed.
        fixed_remaining(75.0)
        assert conversion(*sequences, equivalents=1.5) == pytest.approx(50.0)


class TestRemainingStartingMaterial:

    @pytest.mark.parametrize('measured', [0.0, 42.0, 100.0])
    def test_reports_measured_value_unchanged(self, sequences, fixed_remaining, measured):
        fixed_remaining(measured)
        assert remaining(*sequences) == pytest.approx(measured)

    def test_is_not_clamped_above_one_hundred_percent(self, sequences, fixed_remaining):
        fixed_remaining(150.0)
        assert remaining(*sequences) == pytest.approx(150.0)
