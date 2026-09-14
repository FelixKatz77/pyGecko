'''
Every reader in file_readers returns the centroids exactly as the file holds them, as a Raw_Scans
in the flat layout, plus a metadata dict with the sample name and whatever run-level acquisition
metadata the format carries: start time, polarity and instrument. Nominal-mass binning is not the
readers' job any more (see test_raw_scans).
'''

import base64
import struct
from datetime import datetime, timedelta, timezone

import netCDF4 as nc
import numpy as np
import pytest
from psims.mzml import MzMLWriter

from pygecko.parsers.file_readers import (extract_scans_from_cdf, extract_scans_from_mzml,
                                          extract_scans_from_mzxml)


MASSES = [100.3, 100.4, 201.1, 100.0, 201.4, 201.2]
INTENSITIES = [5.0, 7.0, 3.0, 2.0, 4.0, 6.0]
SCAN_INDEX = [0, 3, 5]


@pytest.fixture
def cdf_file(tmp_path):
    path = tmp_path / 'SAMPLE-A1_MS1Front.cdf'
    times = [10.0, 10.5, 11.0]
    with nc.Dataset(path, 'w') as dataset:
        dataset.test_ionization_polarity = 'Positive Polarity'
        dataset.experiment_date_time_stamp = '20231130185320+0100'
        dataset.createDimension('scan_number', len(times))
        dataset.createDimension('point_number', len(MASSES))
        dataset.createVariable('scan_acquisition_time', 'f8', ('scan_number',))[:] = times
        dataset.createVariable('scan_index', 'i4', ('scan_number',))[:] = SCAN_INDEX
        dataset.createVariable('mass_values', 'f4', ('point_number',))[:] = MASSES
        dataset.createVariable('intensity_values', 'f4', ('point_number',))[:] = INTENSITIES
    return path


def test_cdf_keeps_centroids_unrounded_in_flat_layout(cdf_file):
    raw, metadata = extract_scans_from_cdf(cdf_file)

    np.testing.assert_array_equal(raw.retention_times, [10000.0, 10500.0, 11000.0])
    np.testing.assert_array_equal(raw.scan_index, SCAN_INDEX)
    np.testing.assert_array_equal(raw.mz, np.asarray(MASSES, dtype=np.float32))
    np.testing.assert_array_equal(raw.intensity, INTENSITIES)
    assert metadata == {'SampleName': 'SAMPLE-A1', 'Polarity': 'positive', 'InstrumentName': None,
                        'AcqTime': datetime(2023, 11, 30, 18, 53, 20,
                                            tzinfo=timezone(timedelta(hours=1)))}


def test_cdf_without_acquisition_metadata_leaves_it_none(tmp_path):
    path = tmp_path / 'SAMPLE-B2_MS1Front.cdf'
    with nc.Dataset(path, 'w') as dataset:
        dataset.createDimension('scan_number', 1)
        dataset.createDimension('point_number', 1)
        dataset.createVariable('scan_acquisition_time', 'f8', ('scan_number',))[:] = [1.0]
        dataset.createVariable('scan_index', 'i4', ('scan_number',))[:] = [0]
        dataset.createVariable('mass_values', 'f4', ('point_number',))[:] = [50.0]
        dataset.createVariable('intensity_values', 'f4', ('point_number',))[:] = [1.0]

    _, metadata = extract_scans_from_cdf(path)

    assert metadata == {'SampleName': 'SAMPLE-B2', 'Polarity': None, 'InstrumentName': None,
                        'AcqTime': None}


@pytest.fixture(params=['positive', 'negative'])
def mzml_file(request, tmp_path):
    '''A two-spectrum mzML with a start time, polarity and instrument model, written by psims.'''
    polarity = request.param
    path = tmp_path / 'sample.mzML'
    with MzMLWriter(str(path)) as writer:
        writer.controlled_vocabularies()
        writer.file_description(['MS1 spectrum', 'centroid spectrum'])
        writer.software_list([{'id': 'test', 'version': '0'}])
        writer.instrument_configuration_list([writer.InstrumentConfiguration(
            id='IC1', component_list=[], params=['Agilent instrument model'])])
        writer.data_processing_list([writer.DataProcessing(
            [writer.ProcessingMethod(order=1, software_reference='test',
                                     params=['Conversion to mzML'])], id='DP1')])
        with writer.run(id='SMP-A1', instrument_configuration='IC1',
                        start_time='2023-11-30T18:53:20Z'):
            with writer.spectrum_list(count=2):
                writer.write_spectrum(np.array(MASSES[:3]), np.array(INTENSITIES[:3]),
                                      id='scan=1', polarity=polarity, scan_start_time=1.0,
                                      encoding={'m/z array': np.float64,
                                                'intensity array': np.float32})
                writer.write_spectrum(np.array(MASSES[3:]), np.array(INTENSITIES[3:]),
                                      id='scan=2', polarity=polarity, scan_start_time=1.5,
                                      encoding={'m/z array': np.float64,
                                                'intensity array': np.float32})
    return path, polarity


def test_mzml_keeps_centroids_and_reads_run_metadata(mzml_file):
    path, polarity = mzml_file

    raw, metadata = extract_scans_from_mzml(path)

    np.testing.assert_array_equal(raw.retention_times, [60000.0, 90000.0])
    np.testing.assert_array_equal(raw.scan_index, [0, 3])
    np.testing.assert_array_equal(raw.mz, MASSES)
    np.testing.assert_array_equal(raw.intensity, INTENSITIES)
    assert raw.intensity.dtype == np.float32
    assert metadata == {'SampleName': 'SMP-A1', 'Polarity': polarity,
                        'InstrumentName': 'Agilent instrument model',
                        'AcqTime': datetime(2023, 11, 30, 18, 53, 20, tzinfo=timezone.utc)}


@pytest.fixture
def mzxml_file(tmp_path):
    '''A hand-built two-scan mzXML: network-order float32 m/z-intensity pairs, positive polarity.'''

    def peaks(mzs, intensities):
        pairs = [value for pair in zip(mzs, intensities) for value in pair]
        return base64.b64encode(struct.pack(f'>{len(pairs)}f', *pairs)).decode()

    path = tmp_path / 'SMP-A1.mzXML'
    path.write_text(f'''<?xml version="1.0" encoding="ISO-8859-1"?>
<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">
  <msRun scanCount="2">
    <scan num="1" msLevel="1" peaksCount="3" polarity="+" retentionTime="PT60S">
      <peaks precision="32" byteOrder="network" pairOrder="m/z-int">{peaks(MASSES[:3], INTENSITIES[:3])}</peaks>
    </scan>
    <scan num="2" msLevel="1" peaksCount="3" polarity="+" retentionTime="PT90S">
      <peaks precision="32" byteOrder="network" pairOrder="m/z-int">{peaks(MASSES[3:], INTENSITIES[3:])}</peaks>
    </scan>
  </msRun>
</mzXML>
''')
    return path


def test_mzxml_keeps_centroids_and_reads_polarity(mzxml_file):
    raw, metadata = extract_scans_from_mzxml(mzxml_file)

    np.testing.assert_array_equal(raw.retention_times, [60000.0, 90000.0])
    np.testing.assert_array_equal(raw.scan_index, [0, 3])
    np.testing.assert_array_equal(raw.mz, np.asarray(MASSES, dtype=np.float32))
    np.testing.assert_array_equal(raw.intensity, INTENSITIES)
    assert metadata == {'SampleName': 'SMP-A1', 'Polarity': 'positive', 'InstrumentName': None,
                        'AcqTime': None}
