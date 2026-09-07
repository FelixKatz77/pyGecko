import numpy as np
import pytest
from pygecko.parsers import MS_Base_Parser
from .conftest import fixture_path

# Loads an Agilent .D directory, which only msConvert can turn into mzML.
pytestmark = pytest.mark.msconvert

def test_fid_base_parser():

    test_path = fixture_path('test_injections/FBS-FA-033-A1.D')
    injection = MS_Base_Parser.load_injection(test_path)
    assert injection is not None
    assert injection.sample_name == 'FBS-FA-033-A1'
    assert type(injection.chromatogram) is np.ndarray
    assert injection.chromatogram.shape[0] == 2
    assert injection.detector == 'MS'
    assert injection.peaks is None

if __name__ == '__main__':
    test_fid_base_parser()
