import numpy as np
import pytest
from pygecko.parsers import FID_Base_Parser
from .conftest import fixture_path

def test_fid_base_parser():

    test_path = fixture_path('test_injections/FBS-FA-033-A1.dx_FID1A.CSV')
    injection = FID_Base_Parser.load_injection(test_path, 3.06)
    assert injection is not None
    assert injection.sample_name == 'FBS-FA-033-A1'
    assert injection.solvent_delay == 3.06
    assert type(injection.chromatogram) is np.ndarray
    assert injection.chromatogram.shape[0] == 2
    assert injection.detector == 'FID'
    assert injection.peaks is None
    assert injection.history[0].operation == 'FID_Base_Parser.load_injection'
    assert injection.history[0].parameters['solvent_delay'] == 3.06


def test_integrate_reproduces_the_areas_pick_peaks_computed():

    '''Reconstructing a boarder as index * scan_rate + t0 lands about 1e-15 from the stored axis
    value, which is enough for the searchsorted lookup to pick the neighbouring scan - it did so for
    a quarter of this fixture's 30000 scans, shifting the narrower of its 69 peaks by up to 2.8%.
    Reading boarders straight off the time axis makes the round trip exact.'''

    test_path = fixture_path('test_injections/FBS-FA-033-A1.dx_FID1A.CSV')
    injection = FID_Base_Parser.load_injection(test_path, 2.5)
    injection.pick_peaks()
    expected = {rt: peak.area for rt, peak in injection.peaks.items()}
    assert len(expected) > 1

    injection.integrate()

    for rt, peak in injection.peaks.items():
        assert peak.area == pytest.approx(expected[rt], rel=1e-12)


def test_a_time_range_selects_peaks_inside_the_window():

    '''time_range was converted to scan indices from absolute zero and then applied to a
    chromatogram already truncated at the solvent delay, so this window was really analysed as
    7.5-10.5 min. Only real data has a solvent delay large enough in scans to show it.'''

    test_path = fixture_path('test_injections/FBS-FA-033-A1.dx_FID1A.CSV')
    injection = FID_Base_Parser.load_injection(test_path, 2.5)

    injection.pick_peaks(time_range=(5.0, 8.0))

    assert injection.peaks
    assert 5.0 <= min(injection.peaks) and max(injection.peaks) <= 8.0
    window = injection.processed_chromatogram[0]
    assert 5.0 <= window[0] and window[-1] <= 8.0
