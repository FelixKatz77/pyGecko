'''Injection.flag_peak changes only what it is given.

Called without a flag it appended None to the peak's flags, and called without an analyte it
overwrote the peak's analyte with None, so flagging the internal standard by retention time
(Analysis.fit_calibration_curve) erased the standard's molecule and polyarc then raised.
'''

import pytest

from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.injection.injection import Injection
from pygecko.gc_tools.sequence.gc_sequence import GC_Sequence

from .conftest import make_fid_injection, make_injection, make_peak


class TestFlagPeak:

    def test_without_a_flag_the_flags_are_unchanged(self):
        peak = make_peak(5.0, flags=['standard'])
        make_injection([peak]).flag_peak(5.0)
        assert peak.flags == ['standard']

    def test_without_an_analyte_the_existing_analyte_is_kept(self):
        peak = make_peak(5.0)
        peak.analyte = Analyte(5.0, name='IS')
        make_injection([peak]).flag_peak(5.0)
        assert peak.analyte.name == 'IS'

    def test_a_given_flag_and_analyte_are_assigned(self):
        peak = make_peak(5.0)
        analyte = Analyte(5.0, name='product')
        make_injection([peak]).flag_peak(5.0, flag='checked', analyte=analyte)
        assert peak.flags == ['checked']
        assert peak.analyte is analyte

    def test_flagging_before_peaks_are_picked_raises(self):
        injection = Injection({'SampleName': 'A1'})
        with pytest.raises(ValueError, match='peaks not picked'):
            injection.flag_peak(5.0)


def test_a_sequence_without_picked_peaks_raises_instead_of_warning():
    # GC_Sequence turns a ValueError from a missing standard peak into a warning; an injection
    # whose peaks were never picked is a caller error and must not be reported as a missing peak.
    sequence = GC_Sequence({}, {'A1': make_fid_injection(sample_name='A1')})
    with pytest.raises(ValueError, match='peaks not picked'):
        sequence.set_internal_standard(4.0)
