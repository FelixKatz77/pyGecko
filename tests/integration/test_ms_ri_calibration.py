import numpy as np
from pygecko.parsers import MS_Base_Parser

def test_ms_ri_calibration():
    ms_inject_path = 'fixtures/test_ri_calibration/FKB-FA-060-A1.mzML'
    ms_ri_path = 'fixtures/test_ri_calibration/FKB-FA-060-RI.mzML'

    ms_inject = MS_Base_Parser.load_injection(ms_inject_path, pos=True)
    ri_conf_ms = MS_Base_Parser.load_ri_calibration(ms_ri_path,12, rt=2.255)

    ms_inject.pick_peaks()
    ri_conf_ms.assign_ris(ms_inject)

    alkanes = np.array([
        ('CCCCCCCCCC', 10,  1.558), ('CCCCCCCCCCC', 11,  1.856),
        ('CCCCCCCCCCCC', 12,  2.255), ('CCCCCCCCCCCCC', 13,  2.736),
        ('CCCCCCCCCCCCCC', 14,  3.269), ('CCCCCCCCCCCCCCC', 15,  3.826),
        ('CCCCCCCCCCCCCCCC', 16,  4.388), ('CCCCCCCCCCCCCCCCC', 17,  4.933),
        ('CCCCCCCCCCCCCCCCCC', 18,  5.467), ('CCCCCCCCCCCCCCCCCCC', 19,  5.976),
        ('CCCCCCCCCCCCCCCCCCCC', 20,  6.463), ('CCCCCCCCCCCCCCCCCCCCC', 21,  6.938),
        ('CCCCCCCCCCCCCCCCCCCCCC', 22,  7.389), ('CCCCCCCCCCCCCCCCCCCCCCC', 23,  7.823),
        ('CCCCCCCCCCCCCCCCCCCCCCCC', 24,  8.239), ('CCCCCCCCCCCCCCCCCCCCCCCCC', 25,  8.637),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCC', 26,  9.024), ('CCCCCCCCCCCCCCCCCCCCCCCCCCC', 27,  9.399),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCCCC', 28,  9.757), ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 29, 10.102),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 30, 10.436), ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 31, 10.776),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 32, 11.14 ), ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 33, 11.538),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 34, 11.99 )
    ], dtype=[('smiles', 'U50'), ('c_count', int), ('rt', float)])

    ris = [1201.0, 1271.0, 1730.0, 1893.0]

    assert ms_inject is not None
    assert ri_conf_ms is not None
    assert np.array_equal(ri_conf_ms.alkanes, alkanes)
    assert [round(peak.ri, 0) for peak in ms_inject.peaks.values()] == ris