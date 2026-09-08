'''Loading .pkl files written before the processing history existed.

Pickle is the project's persistence format and .pkl files are coupled to the class layout
(architecture 8), so adding slots to Injection and Analysis_Settings must not strand saved
sequences. Deleting a slot on a live object leaves exactly the state tuple an older file carries, so
that is how the old layout is reproduced here.

The Analysis_Settings case is the one that actually bites: the settings object is nested inside a
pickled injection and restores itself, so Injection.__setstate__ cannot reach it, and without its
own shim the first pop inside Peak_Detection_FID after loading raises AttributeError.
'''

import pickle

import numpy as np

from pygecko.gc_tools.analysis.analysis_settings import Analysis_Settings
from pygecko.gc_tools.injection.fid_injection import FID_Injection

from .conftest import make_injection, make_peak

CHROMATOGRAM = np.array([[0.0, 0.1, 0.2, 0.3], [1.0, 1.0, 1.0, 1.0]])


def round_trip(obj):
    '''Returns the object after a pickle round trip.'''
    return pickle.loads(pickle.dumps(obj))


class TestOldPicklesStillLoad:

    def test_an_injection_pickled_without_a_history_loads_with_an_empty_history(self):
        injection = make_injection([make_peak(5.0)])
        del injection.history
        assert round_trip(injection).history == []

    def test_an_injection_pickled_without_the_nesting_guard_loads_ready_to_record(self):
        injection = make_injection([make_peak(5.0)])
        del injection._recording
        restored = round_trip(injection)
        assert restored._recording is False
        restored.flag_peak(5.0)
        assert len(restored.history) == 1

    def test_settings_pickled_without_resolved_still_accept_pop(self):
        settings = Analysis_Settings(CHROMATOGRAM)
        del settings._resolved
        restored = round_trip(settings)
        assert restored.pop('sn', 1) == 5
        assert restored._resolved == {'sn': 5}

    def test_an_injection_pickled_without_resolved_settings_can_still_be_processed(self):
        injection = FID_Injection({'SampleName': 'SMP-A1'}, CHROMATOGRAM, solvent_delay=0.1)
        del injection.history
        del injection._recording
        del injection.analysis_settings._resolved
        restored = round_trip(injection)
        assert restored.analysis_settings.pop('height', 0) == 0
        assert restored.history == []

    def test_settings_pickled_with_a_stale_indices_range_still_load(self):
        # indices_range was a slot until the window came to be derived from the chromatogram at the
        # call site. Every .pkl written before that carries a value for it, and __setstate__ has to
        # drop names that are no longer slots rather than raise.
        settings = Analysis_Settings(CHROMATOGRAM)
        state = (None, {name: getattr(settings, name) for name in Analysis_Settings.__slots__})
        state[1]['indices_range'] = [0, None]

        restored = Analysis_Settings.__new__(Analysis_Settings)
        restored.__setstate__(state)

        assert not hasattr(restored, 'indices_range')
        assert restored.pop('sn', 1) == 5
        restored.update(time_range=(1.0, 2.0))
        assert restored.time_range == (1.0, 2.0)

    def test_every_other_attribute_survives_the_round_trip(self):
        injection = make_injection([make_peak(5.0)], sample_name='SMP-B7')
        injection.flag_peak(5.0, flag='standard')
        del injection.history
        restored = round_trip(injection)
        assert restored.sample_name == 'SMP-B7'
        assert restored.peaks[5.0].flags == ['standard']

    def test_a_current_injection_round_trips_with_its_history(self):
        injection = make_injection([make_peak(5.0)])
        injection.flag_peak(5.0, flag='standard')
        restored = round_trip(injection)
        assert [step.operation for step in restored.history] == ['Injection.flag_peak']
        assert restored.history[0].parameters['flag'] == 'standard'
