'''
extract_scans_from_cdf reads an ANDI-MS file into the scans DataFrame every MS parser produces: one
row per scan indexed by retention time in milliseconds, one sorted integer column per nominal m/z,
the maximum intensity where a scan has several centroids rounding to the same nominal mass, and 0
where a scan has no centroid at that mass.
'''

import netCDF4 as nc
import numpy as np
import pandas as pd
import pytest

from pygecko.parsers.file_readers import extract_scans_from_cdf


@pytest.fixture
def cdf_file(tmp_path):
    # Scan 1 has 100.3 and 100.4 both rounding to 100 (max wins), scan 2 lacks m/z 100 entirely.
    path = tmp_path / 'SAMPLE-A1_MS1Front.cdf'
    times = [10.0, 10.5, 11.0]
    masses = [100.3, 100.4, 201.1, 100.0, 201.4, 201.2]
    intensities = [5.0, 7.0, 3.0, 2.0, 4.0, 6.0]
    scan_index = [0, 3, 5]
    with nc.Dataset(path, 'w') as dataset:
        dataset.createDimension('scan_number', len(times))
        dataset.createDimension('point_number', len(masses))
        dataset.createVariable('scan_acquisition_time', 'f8', ('scan_number',))[:] = times
        dataset.createVariable('scan_index', 'i4', ('scan_number',))[:] = scan_index
        dataset.createVariable('mass_values', 'f4', ('point_number',))[:] = masses
        dataset.createVariable('intensity_values', 'f4', ('point_number',))[:] = intensities
    return path


def test_scans_are_binned_to_nominal_mass_with_max_and_zero_fill(cdf_file):
    df, sample_name = extract_scans_from_cdf(cdf_file)

    expected = pd.DataFrame([[7.0, 3.0], [2.0, 4.0], [0.0, 6.0]],
                            index=pd.Index([10000.0, 10500.0, 11000.0], name='retention_time'),
                            columns=[100, 201])
    pd.testing.assert_frame_equal(df, expected, check_column_type=False, check_index_type=False)
    assert list(df.columns) == [100, 201]
    assert sample_name == 'SAMPLE-A1'
