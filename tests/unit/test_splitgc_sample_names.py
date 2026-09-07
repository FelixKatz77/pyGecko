'''Tests that SplitGC_Parser warns when the two detectors disagree about sample identity.

The MS side derives sample names from the .cdf filename, the FID side from the acaml SampleName.
When an export makes those diverge, sample_filter matches nothing on one side and the pipeline
silently produces an empty plate. The parser is where the two naming routes meet, so it is where
the divergence is reported.
'''

import pytest

from pygecko.gc_tools.sequence.fid_sequence import FID_Sequence
from pygecko.gc_tools.sequence.ms_sequence import MS_Sequence
from pygecko.parsers.agilent_fid_parser import Agilent_FID_Parser
from pygecko.parsers.ms_base_parser import MS_Base_Parser
from pygecko.parsers.splitgc_parser import SplitGC_Parser


@pytest.fixture
def rslt_dir(tmp_path):
    '''A minimally valid .rslt folder: SplitGC_Parser only checks that AIA/ exists.'''
    path = tmp_path / 'FBS-FB-021-ALL.rslt'
    (path / 'AIA').mkdir(parents=True)
    return path


@pytest.fixture
def stub_parsers(monkeypatch):
    '''Replaces both single-detector parsers so only the coordinator's logic is under test.'''

    def _stub(fid_names, ms_names):
        monkeypatch.setattr(
            Agilent_FID_Parser, 'load_sequence',
            lambda *args, **kwargs: FID_Sequence({}, {name: object() for name in fid_names}))
        monkeypatch.setattr(
            MS_Base_Parser, 'load_sequence',
            lambda *args, **kwargs: MS_Sequence({}, {name: object() for name in ms_names}))

    return _stub


def user_warnings(recwarn):
    return [w for w in recwarn if issubclass(w.category, UserWarning)]


class TestSampleNameDivergence:

    def test_warns_when_no_sample_names_overlap(self, rslt_dir, stub_parsers):
        stub_parsers(['A1', 'A2'], ['20240101-A1', '20240101-A2'])
        with pytest.warns(UserWarning, match='sample names'):
            SplitGC_Parser.load_sequence(str(rslt_dir), 3.0)

    def test_warning_names_both_sets(self, rslt_dir, stub_parsers):
        stub_parsers(['A1'], ['20240101-A1'])
        with pytest.warns(UserWarning) as record:
            SplitGC_Parser.load_sequence(str(rslt_dir), 3.0)
        message = str(record[0].message)
        assert 'A1' in message and '20240101-A1' in message

    def test_still_returns_both_sequences(self, rslt_dir, stub_parsers):
        stub_parsers(['A1'], ['20240101-A1'])
        with pytest.warns(UserWarning):
            fid_sequence, ms_sequence = SplitGC_Parser.load_sequence(str(rslt_dir), 3.0)
        assert list(fid_sequence.injections) == ['A1']
        assert list(ms_sequence.injections) == ['20240101-A1']

    def test_warns_when_one_side_is_empty(self, rslt_dir, stub_parsers):
        '''An empty side is the usual symptom: sample_filter matched nothing on that detector.'''
        stub_parsers(['A1', 'A2'], [])
        with pytest.warns(UserWarning, match='sample names'):
            SplitGC_Parser.load_sequence(str(rslt_dir), 3.0)

    def test_does_not_warn_when_names_overlap(self, rslt_dir, stub_parsers, recwarn):
        stub_parsers(['A1', 'A2'], ['A1', 'A2'])
        SplitGC_Parser.load_sequence(str(rslt_dir), 3.0)
        assert user_warnings(recwarn) == []

    def test_does_not_warn_on_partial_overlap(self, rslt_dir, stub_parsers, recwarn):
        '''Loading one detector's subset of wells stays legal and must not warn.'''
        stub_parsers(['A1', 'A2'], ['A1'])
        SplitGC_Parser.load_sequence(str(rslt_dir), 3.0)
        assert user_warnings(recwarn) == []
