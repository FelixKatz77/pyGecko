'''Tests for internal-standard assignment across a sequence.'''

import pytest

from pygecko.gc_tools.sequence.gc_sequence import GC_Sequence

from .conftest import make_injection, make_peak

IS_RT = 5.00


def build_sequence(with_standard, without_standard=()):
    injections = {}
    for name in with_standard:
        injections[name] = make_injection([make_peak(IS_RT), make_peak(7.0)], sample_name=name)
    for name in without_standard:
        injections[name] = make_injection([make_peak(7.0)], sample_name=name)
    return GC_Sequence({'sequence_name': 'SEQ'}, injections)


class TestSetInternalStandard:

    def test_assigns_standard_to_every_injection(self):
        sequence = build_sequence(['A1', 'A2'])
        sequence.set_internal_standard(IS_RT, name='IS')
        assert all(inj.internal_standard is not None for inj in sequence)
        assert sequence.internal_standard is not None

    def test_warns_instead_of_raising_when_standard_missing(self):
        sequence = build_sequence(['A1'], without_standard=['A2'])
        with pytest.warns(UserWarning, match='No internal standard peak found'):
            sequence.set_internal_standard(IS_RT, name='IS')

    def test_names_the_offending_injections_in_the_warning(self):
        sequence = build_sequence(['A1'], without_standard=['A2', 'A3'])
        with pytest.warns(UserWarning) as record:
            sequence.set_internal_standard(IS_RT, name='IS')
        message = str(record[0].message)
        assert 'A2' in message and 'A3' in message
        assert '2 injection(s)' in message

    def test_remaining_injections_still_get_their_standard(self):
        sequence = build_sequence(['A1'], without_standard=['A2'])
        with pytest.warns(UserWarning):
            sequence.set_internal_standard(IS_RT, name='IS')
        assert sequence['A1'].internal_standard is not None
        assert sequence['A2'].internal_standard is None

    def test_no_warning_when_all_injections_have_the_standard(self, recwarn):
        sequence = build_sequence(['A1', 'A2'])
        sequence.set_internal_standard(IS_RT, name='IS')
        assert not [w for w in recwarn if 'internal standard' in str(w.message)]
