import numpy as np
from pygecko.parsers import FID_Base_Parser, Agilent_FID_Parser


def test_fid_ri_calibration():
    fid_inject_path = 'fixtures/test_injections/FBS-FA-033-A1.dx_FID1A.CSV'
    fid_ri_path = 'fixtures/test_ri_calibration/FID'

    fid_inject = FID_Base_Parser.load_injection(fid_inject_path, 3.0, pos=True)
    ri_conf_fid = Agilent_FID_Parser.load_ri_calibration(fid_ri_path, 2.4, 12, rt=4.570)

    fid_inject.pick_peaks()
    ri_conf_fid.assign_ris(fid_inject)

    alkanes = np.array([
        ('CCCCCCC',  7,  2.426), ('CCCCCCCC',  8,  2.682),
        ('CCCCCCCCC',  9,  3.008), ('CCCCCCCCCC', 10,  3.452),
        ('CCCCCCCCCCC', 11,  3.986), ('CCCCCCCCCCCC', 12,  4.568),
        ('CCCCCCCCCCCCC', 13,  5.164), ('CCCCCCCCCCCCCC', 14,  5.752),
        ('CCCCCCCCCCCCCCC', 15,  6.322), ('CCCCCCCCCCCCCCCC', 16,  6.87 ),
        ('CCCCCCCCCCCCCCCCC', 17,  7.394), ('CCCCCCCCCCCCCCCCCC', 18,  7.896),
        ('CCCCCCCCCCCCCCCCCCC', 19,  8.377), ('CCCCCCCCCCCCCCCCCCCC', 20,  8.837),
        ('CCCCCCCCCCCCCCCCCCCCC', 21,  9.279), ('CCCCCCCCCCCCCCCCCCCCCC', 22,  9.702),
        ('CCCCCCCCCCCCCCCCCCCCCCC', 23, 10.11 ), ('CCCCCCCCCCCCCCCCCCCCCCCC', 24, 10.503),
        ('CCCCCCCCCCCCCCCCCCCCCCCCC', 25, 10.88 ), ('CCCCCCCCCCCCCCCCCCCCCCCCCC', 26, 11.245),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCCC', 27, 11.599), ('CCCCCCCCCCCCCCCCCCCCCCCCCCCC', 28, 11.962),
        ('CCCCCCCCCCCCCCCCCCCCCCCCCCCCC', 29, 12.353)
    ], dtype=[('smiles', 'U50'), ('c_count', int), ('rt', float)])

    ris = [948.0, 998.0, 1048.0, 1076.0, 1172.0, 1207.0, 1262.0, 1290.0, 1314.0, 1335.0, 1353.0, 1372.0, 1391.0, 1411.0,
           1427.0, 1461.0, 1475.0, 1493.0, 1520.0, 1534.0, 1545.0, 1559.0, 1574.0, 1589.0, 1614.0, 1635.0, 1662.0,
           1676.0, 1700.0, 1744.0, 1764.0, 1813.0, 1853.0, 1955.0, 1971.0, 2005.0, 2014.0, 2107.0, 2137.0, 2193.0,
           2267.0, 2333.0, 2369.0, 2418.0, 2484.0, 2497.0, 2550.0, 2566.0, 2585.0, 2608.0, 2626.0, 2636.0, 2680.0,
           2712.0, 2724.0, 2733.0, 2741.0, 2779.0, 2796.0, 2721.0, 2734.0, 2751.0, 2757.0, 2772.0, 2782.0, 2805.0,
           2809.0, 2819.0]

    assert fid_inject is not None
    assert ri_conf_fid is not None
    assert np.array_equal(ri_conf_fid.alkanes, alkanes)
    assert [round(peak.ri, 0) for peak in fid_inject.peaks.values()] == ris

if __name__ == '__main__':
    test_fid_ri_calibration()