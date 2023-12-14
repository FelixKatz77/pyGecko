import numpy as np
from pygecko.parsers import Agilent_FID_Parser
from datetime import datetime

def test_agilent_fid_parser_xy():

    sequence_directory = 'fixtures/test_sequences/FKB-FA-060-RI'
    sequence = Agilent_FID_Parser.load_sequence(sequence_directory, 2.4)
    injection = sequence.injections['FKB-FA-060-RI']
    assert sequence is not None
    assert len(sequence.injections) == 1
    assert sequence.sequence_name == 'FKB-FA-060-RI'
    assert sequence.instrument_name == 'GC-FID'
    assert sequence.internal_standard is None
    assert sequence.detector == 'FID'
    assert injection is not None
    assert injection.sample_name == 'FKB-FA-060-RI'
    assert injection.acq_method == 'Polyarc_110_20_320_2'
    assert injection.acq_time == datetime(2023, 8, 12,  3, 5, 25)
    assert injection.instrument_name == 'GC-FID'
    assert injection.sample_number == '2'
    assert injection.vial_pos == '97'
    assert injection.solvent_delay == 2.4
    assert type(injection.chromatogram) is np.ndarray
    assert injection.chromatogram.shape[0] == 2
    assert injection.detector == 'FID'
    assert injection.peaks is None

def test_agilent_fid_parser_csv():
    sequence_directory = 'fixtures/test_sequences/FBS-FA-034-RI-FID'
    sequence = Agilent_FID_Parser.load_sequence(sequence_directory, 2.7)
    injection = sequence.injections['FBS-FA-034-RI_Std']
    assert sequence is not None
    assert len(sequence.injections) == 1
    assert sequence.sequence_name == 'FBS-FA-034-RI_Std'
    assert sequence.instrument_name == 'GC-FID'
    assert sequence.internal_standard is None
    assert sequence.detector == 'FID'
    assert injection is not None
    assert injection.sample_name == 'FBS-FA-034-RI_Std'
    assert injection.acq_method == 'Polyarc_80_20_310_1'
    assert injection.acq_time == datetime(2023, 12, 11, 18, 4, 55)
    assert injection.instrument_name == 'GC-FID'
    assert injection.sample_number == '2'
    assert injection.vial_pos == '148'
    assert injection.solvent_delay == 2.7
    assert type(injection.chromatogram) is np.ndarray
    assert injection.chromatogram.shape[0] == 2
    assert injection.detector == 'FID'
    assert injection.peaks is None
