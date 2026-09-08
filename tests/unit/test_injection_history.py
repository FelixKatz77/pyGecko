'''Recording of an injection's processing history.

The history exists so a workflow can be reconstructed after the fact - by a reader, or by an agent
building a tool chain from a stored run. That only works if it records what actually happened:
defaulted arguments as the values they took, nested calls once rather than twice, and the thresholds
the algorithms resolved rather than the ones the caller happened to pass. Most real parameters never
appear in a call at all, since Analysis_Settings.pop computes them from the signal, so a history
that only echoed its arguments would be close to useless.
'''

import inspect
import json

import numpy as np
import pytest

import pygecko.gc_tools as gc_tools
import pygecko.parsers as parsers
from pygecko.gc_tools.analysis.analysis_settings import Analysis_Settings
from pygecko.gc_tools.analysis.retention_indices import RI_Calibration
from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.history import Processing_Step
from pygecko.parsers.fid_base_parser import FID_Base_Parser

from .conftest import make_fid_chromatogram, make_fid_injection, make_injection, make_peak

ALKANE_DTYPE = [('smiles', 'U50'), ('c_count', int), ('rt', float)]
IS_RT = 4.0
ANALYTE_RT = 6.0
SOLVENT_DELAY = 1.0


def write_xy(directory, sample_name='SMP-A1'):
    '''Writes a synthetic chromatogram as a tab-delimited .xy file and returns its path.'''
    path = directory / f'{sample_name}.xy'
    np.savetxt(path, np.transpose(make_fid_chromatogram()), delimiter='\t')
    return path


class TestHistoryIsRecorded:

    def test_a_new_injection_starts_with_an_empty_history(self):
        assert make_injection([make_peak(5.0)]).history == []

    def test_a_recorded_call_appends_one_step_naming_the_class_and_method(self):
        injection = make_injection([make_peak(5.0)])
        injection.flag_peak(5.0, flag='standard')
        assert len(injection.history) == 1
        assert injection.history[0].operation == 'Injection.flag_peak'

    def test_defaulted_arguments_are_recorded_with_their_actual_values(self):
        injection = make_injection([make_peak(5.0)])
        injection.flag_peak(5.0)
        assert injection.history[0].parameters == {'rt': 5.0, 'flag': None, 'tolerance': 0.05, 'analyte': None}

    def test_self_is_not_recorded_as_a_parameter(self):
        injection = make_injection([make_peak(5.0)])
        injection.flag_peak(5.0)
        assert 'self' not in injection.history[0].parameters

    def test_keyword_arguments_are_nested_rather_than_flattened(self):
        # Flattening **kwargs into the top level lets a setting silently overwrite a named
        # parameter of the same name, and makes it impossible to tell afterwards where a value came
        # from. Both are defects of the implementation this layer is adapted from.
        injection = make_fid_injection()
        injection.pick_peaks(prominence_fid=500)
        parameters = injection.history[0].parameters
        assert parameters['kwargs'] == {'prominence_fid': 500}
        assert 'prominence_fid' not in parameters

    @pytest.mark.parametrize('switch, value', [('inplace', False), ('inplace', True)])
    def test_boolean_pipeline_switches_are_recorded(self, switch, value):
        injection = make_fid_injection()
        injection.pick_peaks(**{switch: value})
        assert injection.history[0].parameters[switch] is value

    def test_a_switch_on_an_inherited_method_is_recorded(self):
        injection = make_injection([make_peak(5.0)])
        injection.match_rt(5.0, exclude_standard=False)
        assert injection.history[0].parameters['exclude_standard'] is False

    def test_steps_are_appended_in_call_order(self):
        injection = make_injection([make_peak(5.0), make_peak(6.0)])
        injection.flag_peak(5.0)
        injection.match_rt(6.0)
        injection.match_ri(1000)
        assert [step.operation for step in injection.history] == [
            'Injection.flag_peak', 'Injection.match_rt', 'Injection.match_ri']

    def test_a_call_that_raises_records_nothing(self):
        injection = make_injection([make_peak(5.0)])
        with pytest.raises(ValueError):
            injection.set_internal_standard(100.0)
        assert injection.history == []

    def test_an_injection_without_analysis_settings_still_records(self):
        injection = make_injection([make_peak(5.0)])
        assert injection.analysis_settings is None
        injection.flag_peak(5.0)
        assert injection.history[0].resolved == {}


class TestResolvedParametersAreCaptured:

    def test_resolved_reports_the_values_the_algorithms_actually_used(self):
        injection = make_fid_injection()
        injection.pick_peaks()
        resolved = injection.history[0].resolved
        assert set(resolved) >= {'indices_range', 'prominence_fid', 'width', 'height',
                                 'boarder_threshold', 'boarder_window'}
        # Computed from the signal, never passed by the caller, and invisible without this capture.
        assert resolved['prominence_fid'] > 0

    def test_a_configured_setting_appears_in_resolved_instead_of_its_default(self):
        injection = make_fid_injection()
        injection.pick_peaks(prominence_fid=42)
        assert injection.history[0].resolved['prominence_fid'] == 42

    def test_resolved_is_cleared_between_steps(self):
        injection = make_fid_injection()
        injection.baseline_correction()
        injection.match_rt(ANALYTE_RT) if injection.peaks else None
        injection.pick_peaks()
        assert 'prominence_fid' not in injection.history[0].resolved
        assert 'prominence_fid' in injection.history[-1].resolved

    def test_resolved_is_not_a_settable_setting(self):
        settings = Analysis_Settings(make_fid_chromatogram())
        with pytest.raises(KeyError):
            settings.update(_resolved={})

    def test_resolved_is_not_a_poppable_setting(self):
        settings = Analysis_Settings(make_fid_chromatogram())
        with pytest.raises(KeyError):
            settings.pop('_resolved', 1)

    def test_a_falsy_configured_setting_falls_through_to_the_default_and_the_default_is_recorded(self):
        # Analysis_Settings.pop tests truthiness, so a legitimately configured 0 is treated as
        # unset (architecture 5, a known deviation). Pinned rather than fixed: resolved must report
        # the value the caller actually received, whatever that value is.
        settings = Analysis_Settings(make_fid_chromatogram())
        settings.update(sn=0)
        assert settings.pop('sn', 5) == 5
        assert settings._resolved['sn'] == 5


class TestNestedCallsRecordOnce:

    def test_pick_peaks_does_not_record_its_internal_baseline_correction(self):
        injection = make_fid_injection()
        injection.pick_peaks()
        assert [step.operation for step in injection.history] == ['FID_Injection.pick_peaks']

    def test_the_nested_baseline_correction_settings_appear_in_the_outer_step(self):
        injection = make_fid_injection()
        injection.pick_peaks()
        assert set(injection.history[0].resolved) >= {'savgol_window', 'max_half_window'}

    def test_set_internal_standard_does_not_record_its_internal_flag_peak(self):
        injection = make_injection([make_peak(5.0)])
        injection.set_internal_standard(5.0, name='dodecane')
        assert [step.operation for step in injection.history] == ['Injection.set_internal_standard']

    def test_an_explicitly_called_baseline_correction_is_recorded(self):
        # The guard suppresses nesting, not repetition.
        injection = make_fid_injection()
        injection.baseline_correction()
        injection.pick_peaks()
        assert [step.operation for step in injection.history] == [
            'FID_Injection.baseline_correction', 'FID_Injection.pick_peaks']

    def test_the_guard_is_released_after_a_failing_call(self):
        injection = make_injection([make_peak(5.0)])
        with pytest.raises(ValueError):
            injection.set_internal_standard(100.0)
        injection.flag_peak(5.0)
        assert [step.operation for step in injection.history] == ['Injection.flag_peak']


class TestLoadProvenance:

    def test_record_step_appends_an_encoded_step(self):
        injection = make_injection([make_peak(5.0)])
        injection.record_step('Some_Parser.load_injection', {'path': 'a.xy', 'time_range': (1, 2)})
        assert injection.history[0].operation == 'Some_Parser.load_injection'
        assert injection.history[0].parameters == {'path': 'a.xy', 'time_range': [1, 2]}

    def test_the_parser_records_the_load_as_the_first_step(self, tmp_path):
        injection = FID_Base_Parser.load_injection(write_xy(tmp_path), SOLVENT_DELAY)
        assert injection.history[0].operation == 'FID_Base_Parser.load_injection'
        assert injection.history[0].parameters['solvent_delay'] == SOLVENT_DELAY

    def test_the_load_step_names_the_source_file(self, tmp_path):
        path = write_xy(tmp_path)
        injection = FID_Base_Parser.load_injection(path, SOLVENT_DELAY)
        assert injection.history[0].parameters['xy_file'] == str(path)

    def test_the_load_step_precedes_the_processing_steps(self, tmp_path):
        injection = FID_Base_Parser.load_injection(write_xy(tmp_path), SOLVENT_DELAY)
        injection.pick_peaks()
        assert [step.operation for step in injection.history] == [
            'FID_Base_Parser.load_injection', 'FID_Injection.pick_peaks']


class TestRiAssignmentIsRecordedOnTheInjection:

    @pytest.fixture
    def calibration(self):
        '''An RI_Calibration built from a fitted ladder, bypassing peak detection on a raw file.'''
        calibration = object.__new__(RI_Calibration)
        calibration.alkanes = np.array(
            [('CCCCCCCCCC', 10, 4.0), ('CCCCCCCCCCC', 11, 5.0),
             ('CCCCCCCCCCCC', 12, 6.0), ('CCCCCCCCCCCCC', 13, 7.0)], dtype=ALKANE_DTYPE)
        calibration.gradient = 100.0
        calibration.intercept = 600.0
        return calibration

    def test_assign_ris_records_a_step_on_the_injection_it_mutates(self, calibration):
        injection = make_injection([make_peak(5.5)])
        calibration.assign_ris(injection)
        assert injection[5.5].ri is not None
        assert [step.operation for step in injection.history] == ['RI_Calibration.assign_ris']

    def test_the_assign_ris_step_carries_the_calibration_needed_to_recompute_the_indices(self, calibration):
        injection = make_injection([make_peak(5.5)])
        calibration.assign_ris(injection)
        parameters = injection.history[0].parameters
        assert parameters['alignment'] is False
        assert parameters['gradient'] == 100.0
        assert parameters['intercept'] == 600.0
        assert parameters['alkanes'][0] == {'smiles': 'CCCCCCCCCC', 'c_count': 10, 'rt': 4.0}


class TestHistoryExport:

    def test_history_to_json_returns_every_step_with_its_operation_and_parameters(self):
        injection = make_injection([make_peak(5.0), make_peak(6.0)])
        injection.flag_peak(5.0, flag='standard')
        injection.match_rt(6.0)
        document = json.loads(injection.history_to_json())
        assert [step['operation'] for step in document] == ['Injection.flag_peak', 'Injection.match_rt']
        assert document[0]['parameters']['flag'] == 'standard'

    def test_history_to_json_writes_the_document_when_a_path_is_given(self, tmp_path):
        injection = make_injection([make_peak(5.0)])
        injection.flag_peak(5.0)
        path = tmp_path / 'history.json'
        injection.history_to_json(path=str(path))
        assert json.loads(path.read_text())[0]['operation'] == 'Injection.flag_peak'


class TestHistoryIsSufficientForReplay:

    @pytest.fixture
    def processed(self, tmp_path):
        '''An injection carrying a full FID workflow, from the raw file to a quantified yield.'''
        injection = FID_Base_Parser.load_injection(write_xy(tmp_path), SOLVENT_DELAY)
        injection.pick_peaks(prominence_fid=100)
        injection.set_internal_standard(IS_RT, name='dodecane', smiles='CCCCCCCCCCCC')
        # No integrate() here: it re-integrates the raw chromatogram, indexing the minute-valued
        # boarders pick_peaks produced as if they were scan indices, and raises on real data. That
        # is a pre-existing defect, unrelated to the history, and left alone deliberately.
        injection.flag_peak(ANALYTE_RT, analyte=Analyte(ANALYTE_RT, smiles='c1ccccc1'))
        injection.quantify(min(injection.peaks, key=lambda rt: abs(rt - ANALYTE_RT)))
        return injection

    @staticmethod
    def resolve(operation):
        '''Returns the callable a recorded operation names, looked up as gc_tools then parsers.'''
        class_name, method_name = operation.split('.')
        for module in (gc_tools, parsers):
            if hasattr(module, class_name):
                return getattr(getattr(module, class_name), method_name)
        raise AssertionError(f'{operation} does not resolve to a class in gc_tools or parsers.')

    def test_every_recorded_operation_resolves_to_a_callable(self, processed):
        assert processed.history
        for step in processed.history:
            assert callable(self.resolve(step.operation)), step.operation

    def test_every_recorded_parameter_set_is_accepted_by_its_operations_signature(self, processed):
        for step in processed.history:
            function = self.resolve(step.operation)
            signature = inspect.signature(function)
            parameters = dict(step.parameters)
            arguments = {**parameters.pop('kwargs', {}), **parameters}
            if 'self' in signature.parameters:
                signature.bind(processed, **arguments)
            else:
                signature.bind(**arguments)

    def test_a_full_workflow_history_round_trips_through_json_with_nothing_unserializable(self, processed):
        document = json.loads(processed.history_to_json())
        assert [step['operation'] for step in document] == [
            'FID_Base_Parser.load_injection', 'FID_Injection.pick_peaks',
            'Injection.set_internal_standard', 'Injection.flag_peak',
            'FID_Injection.quantify']
        assert 'unserializable' not in json.dumps(document)
