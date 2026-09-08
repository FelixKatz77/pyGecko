import functools
import inspect
from datetime import datetime, timezone

import numpy as np

from pygecko import __version__
from pygecko.gc_tools.analyte import Analyte


class Processing_Step:

    '''
    A single operation applied to an injection, together with the parameters it ran with.

    Attributes:
        operation (str): Qualified name of the operation as Class.method, e.g.
            'FID_Injection.pick_peaks'. This is the pair a replay driver dispatches on.
        parameters (dict): Encoded arguments the operation was called with, defaults applied.
            Variable keyword arguments stay nested under their own parameter name rather than being
            flattened into the top level, so a setting can never collide with a named parameter.
        resolved (dict): Encoded values Analysis_Settings.pop handed to the algorithms during the
            call, keyed by setting name. Empty for operations that read no settings.
        timestamp (str): ISO 8601 UTC time the step was recorded.
        version (str): pyGecko version that recorded the step.
    '''

    operation: str
    parameters: dict
    resolved: dict
    timestamp: str
    version: str

    __slots__ = 'operation', 'parameters', 'resolved', 'timestamp', 'version'

    def __init__(self, operation:str, parameters:dict, resolved:dict|None=None):
        self.operation = operation
        self.parameters = Processing_Step.__encode(parameters)
        self.resolved = Processing_Step.__encode(resolved or {})
        self.timestamp = datetime.now(timezone.utc).isoformat()
        self.version = __version__

    def __str__(self) -> str:

        '''
        Returns the step as the call it recorded.
        '''

        arguments = ', '.join(f'{name}={value!r}' for name, value in self.parameters.items())
        return f'{self.operation}({arguments})'

    def to_dict(self) -> dict:

        '''
        Returns the step as a JSON-serializable dictionary.

        Returns:
            dict: The step's operation, parameters, resolved settings, timestamp and version.
        '''

        return {'operation': self.operation, 'parameters': self.parameters, 'resolved': self.resolved,
                'timestamp': self.timestamp, 'version': self.version}

    @staticmethod
    def __encode(value):

        '''
        Returns a JSON-serializable representation of a recorded value.

        Tuples become lists, since JSON has neither tuples nor sets; time_range is the one setting
        Analysis_Settings type-checks as a tuple, so a replay driver has to restore that type itself.
        Structured arrays become lists of records keyed by field name, which is what lets the RI
        calibration's alkane ladder survive with its columns named. A value with no JSON
        representation - a callable above all - is encoded as an explicit marker naming its type and
        repr, never dropped: a missing parameter is indistinguishable from one that was never
        passed, which is exactly what makes a silently filtered history unreplayable.

        Args:
            value: Value to encode.

        Returns:
            The JSON-serializable representation of the value.
        '''

        # numpy scalars first: np.float64 is a subclass of float, so a scalar check would catch it
        # and leave a numpy type in the record that json cannot serialize.
        if isinstance(value, np.generic):
            return value.item()
        if value is None or isinstance(value, (bool, int, float, str)):
            return value
        if isinstance(value, dict):
            return {str(key): Processing_Step.__encode(item) for key, item in value.items()}
        if isinstance(value, (list, tuple, set, frozenset)):
            return [Processing_Step.__encode(item) for item in value]
        if isinstance(value, np.ndarray):
            if value.dtype.names:
                return [dict(zip(value.dtype.names, Processing_Step.__encode(record)))
                        for record in value.tolist()]
            return Processing_Step.__encode(value.tolist())
        if isinstance(value, Analyte):
            return {'rt': value.rt, 'name': value.name, 'smiles': value.smiles}
        return {'unserializable': type(value).__name__, 'repr': repr(value)}


def records_processing(func):

    '''
    Decorator appending a Processing_Step to an injection's history for each call.

    The step records the call's arguments with defaults applied, plus every setting the algorithms
    resolved out of the injection's Analysis_Settings while it ran. Only the outermost decorated
    call on an injection records: FID_Injection.pick_peaks calls baseline_correction when no
    processed chromatogram exists yet, and set_internal_standard calls flag_peak, so recording the
    inner call would put a step in the history the caller never asked for and a replay driver would
    then execute twice. The guard suppresses nesting, not repetition - an explicit
    baseline_correction followed by pick_peaks still records two steps.

    Nothing is recorded when the call raises, since the state the step describes was never reached.

    Args:
        func (callable): Injection method to record.

    Returns:
        callable: The wrapped method.
    '''

    # Bound once here rather than per call: flag_peak runs once per peak in
    # RI_Calibration.__identify_alkanes and once per well in Analysis.
    signature = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self._recording:
            return func(self, *args, **kwargs)
        settings = self.analysis_settings
        # Cleared once, at the outer call only. This is what lets the single pick_peaks step report
        # both the nested baseline-correction settings and its own; clearing on every call would
        # wipe the former before the outer call snapshots them.
        if settings is not None:
            settings._resolved.clear()
        self._recording = True
        try:
            result = func(self, *args, **kwargs)
        finally:
            self._recording = False
        resolved = settings._resolved if settings is not None else {}
        parameters = _bind_parameters(signature, (self,) + args, kwargs)
        self.history.append(Processing_Step(func.__qualname__, parameters, resolved))
        return result

    return wrapper


def _bind_parameters(signature:inspect.Signature, args:tuple, kwargs:dict) -> dict:

    '''
    Returns a call's arguments keyed by parameter name, with defaults applied and self excluded.

    Variable keyword arguments are left nested under their own parameter name, so a setting passed
    through **kwargs can never silently overwrite a named parameter of the same name.

    Args:
        signature (inspect.Signature): Signature of the recorded method.
        args (tuple): Positional arguments of the call, including self.
        kwargs (dict): Keyword arguments of the call.

    Returns:
        dict: The bound arguments, excluding self.
    '''

    bound = signature.bind(*args, **kwargs)
    bound.apply_defaults()
    parameters = dict(bound.arguments)
    parameters.pop('self', None)
    return parameters
