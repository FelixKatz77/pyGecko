'''
MS_Base_Parser builds an MS_Injection from a real OpenChrom-exported mzML: the centroids as read
next to the nominal-mass matrix derived from them, and the run metadata the file carries.
'''

from datetime import datetime, timedelta, timezone

import numpy as np

from pygecko.parsers import MS_Base_Parser

from .conftest import fixture_path


REAL_MZML = fixture_path('test_ri_calibration', 'FKB-FA-060-A1.mzML')


def test_load_injection_keeps_raw_scans_next_to_the_nominal_matrix():
    injection = MS_Base_Parser.load_injection(REAL_MZML)

    assert len(injection.raw_scans.retention_times) == 1878 == len(injection.scans)
    # The vendor file reports 0.05-Da centroids; the raw arrays keep them unrounded.
    assert not np.all(injection.raw_scans.mz == np.round(injection.raw_scans.mz))
    np.testing.assert_array_equal(injection.scans.index, injection.raw_scans.retention_times)


def test_load_injection_tic_equals_the_sum_of_the_centroids():
    injection = MS_Base_Parser.load_injection(REAL_MZML)

    expected = [intensity.sum() for _, _, intensity in injection.raw_scans.spectra()]
    np.testing.assert_allclose(injection.chromatogram[1], expected, rtol=1e-6)


def test_load_injection_reads_the_run_metadata():
    injection = MS_Base_Parser.load_injection(REAL_MZML)

    assert injection.sample_name == 'FKB-FA-060-A1'
    assert injection.acq_time == datetime(2023, 8, 10, 2, 27, tzinfo=timezone(timedelta(hours=2)))
    assert injection.polarity is None  # the OpenChrom export carries no polarity term
    assert injection.instrument_name is None  # its instrument configuration is empty
