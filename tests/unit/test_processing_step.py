'''Encoding behaviour of a single recorded processing step.

A history is only worth keeping if it is faithful. The reference implementation this layer is
adapted from filtered recorded arguments through a scalar whitelist and silently dropped anything
else, so a step could describe a call it did not make: a dropped parameter is indistinguishable from
one that was never passed. These tests pin the opposite contract - every argument is either encoded
losslessly or recorded as an explicit marker naming what could not be encoded.
'''

import json

import numpy as np
import pytest

from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.history import Processing_Step

ALKANE_DTYPE = [('smiles', 'U50'), ('c_count', int), ('rt', float)]


class TestProcessingStep:

    def test_a_step_records_its_operation_parameters_and_resolved_settings(self):
        step = Processing_Step('FID_Injection.pick_peaks', {'inplace': True}, {'prominence_fid': 12.5})
        assert step.operation == 'FID_Injection.pick_peaks'
        assert step.parameters == {'inplace': True}
        assert step.resolved == {'prominence_fid': 12.5}

    def test_a_step_without_resolved_settings_records_an_empty_mapping(self):
        step = Processing_Step('Injection.flag_peak', {'rt': 5.0})
        assert step.resolved == {}

    def test_a_step_records_when_it_ran_and_which_version_recorded_it(self):
        from pygecko import __version__

        step = Processing_Step('Injection.flag_peak', {'rt': 5.0})
        assert step.version == __version__
        assert step.timestamp.startswith(str(np.datetime64('now', 'Y')))

    def test_a_step_serializes_to_a_json_document(self):
        step = Processing_Step('FID_Injection.pick_peaks', {'kwargs': {'width': 3}}, {'height': 0})
        document = json.loads(json.dumps(step.to_dict()))
        assert document['operation'] == 'FID_Injection.pick_peaks'
        assert document['parameters'] == {'kwargs': {'width': 3}}
        assert document['resolved'] == {'height': 0}
        assert set(document) == {'operation', 'parameters', 'resolved', 'timestamp', 'version'}

    def test_a_step_reads_as_the_call_it_recorded(self):
        step = Processing_Step('Injection.flag_peak', {'rt': 5.0, 'flag': 'standard'})
        assert str(step) == "Injection.flag_peak(rt=5.0, flag='standard')"


class TestParameterEncoding:

    @pytest.mark.parametrize('value', [None, True, False, 0, 3, 3.5, 'a'])
    def test_json_scalars_pass_through_unchanged(self, value):
        assert Processing_Step('op', {'x': value}).parameters['x'] == value

    def test_numpy_scalars_become_python_scalars(self):
        encoded = Processing_Step('op', {'x': np.mean(np.array([1.0, 2.0]))}).parameters['x']
        assert encoded == 1.5
        assert type(encoded) is float

    def test_numpy_integers_become_python_integers(self):
        encoded = Processing_Step('op', {'x': np.int64(7)}).parameters['x']
        assert encoded == 7
        assert type(encoded) is int

    def test_tuples_are_encoded_as_lists(self):
        # JSON has no tuple. time_range is the one setting Analysis_Settings type-checks as a
        # tuple, so a replay driver has to restore that type itself.
        assert Processing_Step('op', {'time_range': (1.0, 9.0)}).parameters['time_range'] == [1.0, 9.0]

    def test_nested_containers_are_encoded_recursively(self):
        parameters = {'kwargs': {'time_range': (1, 2), 'widths': [np.float64(0.5)]}}
        assert Processing_Step('op', parameters).parameters == {'kwargs': {'time_range': [1, 2], 'widths': [0.5]}}

    def test_analytes_are_encoded_as_rt_name_and_smiles(self):
        analyte = Analyte(5.0, name='dodecane', smiles='CCCCCCCCCCCC')
        assert Processing_Step('op', {'analyte': analyte}).parameters['analyte'] == {
            'rt': 5.0, 'name': 'dodecane', 'smiles': 'CCCCCCCCCCCC'}

    def test_structured_arrays_are_encoded_as_records_keyed_by_field_name(self):
        alkanes = np.array([('CCCCCCCCCC', 10, 4.5), ('CCCCCCCCCCC', 11, 5.5)], dtype=ALKANE_DTYPE)
        assert Processing_Step('op', {'alkanes': alkanes}).parameters['alkanes'] == [
            {'smiles': 'CCCCCCCCCC', 'c_count': 10, 'rt': 4.5},
            {'smiles': 'CCCCCCCCCCC', 'c_count': 11, 'rt': 5.5},
        ]

    def test_plain_arrays_are_encoded_as_nested_lists(self):
        assert Processing_Step('op', {'boarders': np.array([[1.0, 2.0], [3.0, 4.0]])}).parameters['boarders'] \
               == [[1.0, 2.0], [3.0, 4.0]]

    def test_callables_are_recorded_as_an_explicit_marker(self):
        encoded = Processing_Step('op', {'func': lambda x: x + 1}).parameters['func']
        assert encoded['unserializable'] == 'function'
        assert 'lambda' in encoded['repr']

    def test_arbitrary_objects_are_recorded_as_an_explicit_marker(self):
        encoded = Processing_Step('op', {'thing': object()}).parameters['thing']
        assert encoded['unserializable'] == 'object'
        assert encoded['repr'].startswith('<object object')

    def test_a_marker_is_still_serializable(self):
        step = Processing_Step('op', {'func': lambda x: x})
        json.dumps(step.to_dict())

    def test_parameters_are_snapshotted_at_record_time(self):
        # Encoding eagerly is what keeps a step honest: holding the live Analyte would let a later
        # assignment rewrite a record of something that already happened.
        analyte = Analyte(5.0, name='dodecane')
        step = Processing_Step('op', {'analyte': analyte})
        analyte.name = 'undecane'
        assert step.parameters['analyte']['name'] == 'dodecane'

    def test_resolved_settings_are_encoded_too(self):
        assert Processing_Step('op', {}, {'prominence_fid': np.float64(1284.7)}).resolved == {'prominence_fid': 1284.7}
