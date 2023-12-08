import numpy as np
from parsers import FID_Base_Parser, Agilent_FID_Parser

def test_fid_base_parser():

    test_path = 'fixtures/test_injections/FBS-FA-033-A1.dx_FID1A.CSV'
    injection = FID_Base_Parser.load_injection(test_path, 3.06)
    assert injection is not None
    assert injection.sample_name == 'FBS-FA-033-A1'
    assert injection.solvent_delay == 3.06
    assert type(injection.chromatogram) is np.ndarray
    assert injection.chromatogram.shape[0] == 2
    assert injection.detector == 'FID'
    assert injection.peaks is None




