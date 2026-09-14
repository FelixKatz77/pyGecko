'''
Raw_Scans holds the centroids of an injection exactly as the reader delivered them, in the flat
ANDI-style layout (one m/z and intensity array for the whole run, a start offset per scan), and is
the one place they are binned to the nominal-mass matrix every MS parser produces.
'''

import numpy as np
import pandas as pd
import pytest

from pygecko.gc_tools.injection.raw_scans import Raw_Scans


@pytest.fixture
def raw():
    # Scan 1: 100.3 and 100.4 both round to 100. Scan 2 is empty. Scan 3 has one centroid.
    return Raw_Scans(
        retention_times=np.array([10000.0, 10500.0, 11000.0]),
        scan_index=np.array([0, 3, 3]),
        mz=np.array([100.3, 100.4, 201.1, 201.4]),
        intensity=np.array([5.0, 7.0, 3.0, 6.0], dtype=np.float32),
    )


def test_to_nominal_matrix_sums_split_centroids_and_zero_fills(raw):
    expected = pd.DataFrame([[12.0, 3.0], [0.0, 0.0], [0.0, 6.0]],
                            index=pd.Index([10000.0, 10500.0, 11000.0], name='retention_time'),
                            columns=[100, 201])

    pd.testing.assert_frame_equal(raw.to_nominal_matrix(), expected,
                                  check_column_type=False, check_index_type=False)


def test_to_nominal_matrix_has_sorted_integer_columns_and_float_values():
    raw = Raw_Scans(np.array([1.0]), np.array([0]), np.array([300.2, 41.0, 41.4]),
                    np.array([1.0, 2.0, 3.0], dtype=np.float32))

    matrix = raw.to_nominal_matrix()

    assert list(matrix.columns) == [41, 300]
    assert matrix.dtypes.map(lambda dtype: dtype == np.float64).all()


def test_spectra_yields_one_slice_per_scan_in_the_original_order(raw):
    spectra = list(raw.spectra())

    assert [rt for rt, _, _ in spectra] == [10000.0, 10500.0, 11000.0]
    np.testing.assert_array_equal(spectra[0][1], [100.3, 100.4, 201.1])
    np.testing.assert_array_equal(spectra[0][2], [5.0, 7.0, 3.0])
    assert spectra[1][1].size == 0
    np.testing.assert_array_equal(spectra[2][1], [201.4])
    assert spectra[2][2].dtype == np.float32
