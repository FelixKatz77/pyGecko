import numpy as np
from pygecko.parsers import Agilent_MS_Parser
from datetime import datetime

def test_agilent_fid_parser():

    sequence_directory = 'fixtures/test_sequences/FKB-FA-059-II'
    sequence = Agilent_MS_Parser.load_sequence(sequence_directory)
    injection = sequence.injections['FKB-FA-059-II-C12']
    assert sequence is not None
    assert len(sequence.injections) == 3
    assert sequence.sequence_name == 'FKB-FA-059-II'
    assert sequence.instrument_name == 'GCMS 4'
    assert sequence.internal_standard is None
    assert sequence.detector == 'MS'
    assert injection is not None
    assert injection.sample_name == 'FKB-FA-059-II-C12'
    assert injection.acq_method == '110_20_320_2_1,4_MINSD'
    assert injection.instrument_name == 'GCMS 4'
    assert injection.vial_pos == '12'
    assert type(injection.chromatogram) is np.ndarray
    assert injection.chromatogram.shape[0] == 2
    assert injection.detector == 'MS'
    assert injection.peaks is None

if __name__ == '__main__':
    test_agilent_fid_parser()