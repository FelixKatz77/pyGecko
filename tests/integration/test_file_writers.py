'''Round-trip tests for the export writers.

Every test here writes a file and reads it back, with the reader pyGecko already ships or with an
independent parser, so the success criterion is what a reader sees rather than a hand-checked byte
layout. An MS export carries the centroids as read, so fidelity is checked against the original
vendor mzML as well as against pyGecko's own reader.
'''

from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET

import netCDF4 as nc
import numpy as np
import pandas as pd
import pytest
from pyteomics import mzml

from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.injection.ms_injection import MS_Injection
from pygecko.gc_tools.sequence.fid_sequence import FID_Sequence
from pygecko.gc_tools.sequence.ms_sequence import MS_Sequence
from pygecko.parsers.agilent_fid_parser import Agilent_FID_Parser
from pygecko.parsers.file_readers import extract_scans_from_mzml
from pygecko.parsers.file_writers import (write_injection_to_cdf, write_injection_to_mzml,
                                          write_sequence_to_cdf, write_sequence_to_mzml)
from pygecko.parsers.ms_base_parser import MS_Base_Parser

from .conftest import fixture_path
from ..unit.conftest import make_fid_chromatogram, make_scans


REAL_MZML = fixture_path('test_ri_calibration', 'FKB-FA-060-A1.mzML')


def read_scans(path):
    '''Reads an mzML back through pyGecko's reader into the nominal matrix and sample name.'''
    raw, metadata = extract_scans_from_mzml(path)
    return raw.to_nominal_matrix(), metadata['SampleName']


def read_spectra(path):
    '''Reads every spectrum of an mzML with pyteomics, independently of pymzml and psims.'''
    with mzml.read(str(path)) as reader:
        return list(reader)


@pytest.fixture
def ms_injection():
    '''An MS_Injection holding a nominal matrix only, as a file saved before raw scans existed.'''

    scans = make_scans(
        [89575.0, 89926.0, 90278.0],
        [{40: 1000.0, 41: 250.5}, {40: 900.0, 78: 12.0}, {41: 5.0, 78: 3000.25}],
    )
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    return MS_Injection({'SampleName': 'SMP-A1'}, chromatogram, None, scans)


@pytest.fixture
def real_injection():
    return MS_Base_Parser.load_injection(REAL_MZML)


def test_mzml_export_matches_the_source_file_spectrum_for_spectrum(real_injection, tmp_path):
    '''What the vendor mzML holds is what the export holds: unrounded m/z, intensities, times.'''

    out = tmp_path / 'FKB-FA-060-A1.mzML'
    write_injection_to_mzml(real_injection, out)

    source, export = read_spectra(REAL_MZML), read_spectra(out)

    assert len(export) == len(source) == 1878
    for original, written in zip(source, export):
        np.testing.assert_array_equal(written['m/z array'], original['m/z array'])
        np.testing.assert_array_equal(written['intensity array'], original['intensity array'])
        assert written['intensity array'].dtype == original['intensity array'].dtype
        assert (written['scanList']['scan'][0]['scan start time']
                == original['scanList']['scan'][0]['scan start time'])
        assert written['total ion current'] == pytest.approx(original['total ion current'])


def test_mzml_export_writes_only_the_metadata_the_source_had(real_injection, tmp_path):
    '''The OpenChrom fixture has a start time but no polarity: neither is invented nor dropped.'''

    out = tmp_path / 'FKB-FA-060-A1.mzML'
    write_injection_to_mzml(real_injection, out)

    run = ET.parse(out).getroot().find('.//{*}run')
    assert datetime.fromisoformat(run.attrib['startTimeStamp']) == real_injection.acq_time
    spectrum = read_spectra(out)[0]
    assert 'positive scan' not in spectrum and 'negative scan' not in spectrum


def test_mzml_export_reads_back_through_the_parser_unchanged(real_injection, tmp_path):
    '''Raw scans, matrix and metadata survive pyGecko's own read of the export.'''

    out = tmp_path / 'FKB-FA-060-A1.mzML'
    write_injection_to_mzml(real_injection, out)

    round_tripped = MS_Base_Parser.load_injection(out)

    np.testing.assert_array_equal(round_tripped.raw_scans.mz, real_injection.raw_scans.mz)
    np.testing.assert_array_equal(round_tripped.raw_scans.intensity,
                                  real_injection.raw_scans.intensity)
    pd.testing.assert_frame_equal(round_tripped.scans, real_injection.scans)
    assert round_tripped.sample_name == 'FKB-FA-060-A1'
    assert round_tripped.acq_time == real_injection.acq_time


def test_mzml_writes_polarity_and_instrument_when_known(tmp_path):
    '''Metadata a vendor conversion provides (msConvert on a .D) is written as mzML terms.'''

    scans = make_scans([1000.0, 1500.0, 2000.0], [{40: 1.0}, {40: 2.0}, {40: 3.0}])
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    injection = MS_Injection({'SampleName': 'SMP-A1', 'InstrumentName': 'GCMS 4',
                              'Polarity': 'positive',
                              'AcqTime': datetime(2023, 11, 30, 18, 53, 20, tzinfo=timezone.utc)},
                             chromatogram, None, scans)

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(injection, out)

    assert 'positive scan' in read_spectra(out)[0]
    _, metadata = extract_scans_from_mzml(out)
    assert metadata['InstrumentName'] == 'GCMS 4'
    assert metadata['Polarity'] == 'positive'
    assert metadata['AcqTime'] == injection.acq_time


def test_mzml_round_trip_preserves_scans(ms_injection, tmp_path):
    '''An injection without raw scans is written from its matrix and reads back exactly.'''

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    scans, sample_name = read_scans(out)

    pd.testing.assert_frame_equal(scans, ms_injection.scans)
    assert sample_name == 'SMP-A1'


def test_mzml_written_from_the_matrix_declares_the_binning(ms_injection, tmp_path):
    '''Without raw scans the file holds nominal masses, and its processing history says so.'''

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    names = [param.attrib['name'] for param in
             ET.parse(out).getroot().findall('.//{*}processingMethod/{*}userParam')]
    assert names == ['nominal mass binning']


def test_mzml_is_readable_by_an_independent_parser(ms_injection, tmp_path):
    '''pyteomics agrees with pymzml on the written file, so the output is not pymzml-specific.'''

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    spectra = read_spectra(out)

    assert len(spectra) == len(ms_injection.scans)
    assert spectra[0]['ms level'] == 1
    # Zeros are dropped on write and restored by the reader's zero fill, so the written arrays
    # only carry the m/z actually present in that scan.
    assert list(spectra[0]['m/z array']) == [40.0, 41.0]


def test_mzml_drops_zero_intensities(ms_injection, tmp_path):
    '''The dense zero-filled matrix is written sparsely.'''

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    lengths = [len(s['m/z array']) for s in read_spectra(out)]

    assert lengths == [2, 2, 2]
    assert ms_injection.scans.shape[1] == 3


def test_mzml_sanitises_a_sample_name_that_is_not_a_valid_xml_id(tmp_path):
    '''A sample name with spaces or a leading digit still produces a readable file.'''

    scans = make_scans([1000.0, 1500.0, 2000.0], [{40: 1.0}, {40: 2.0}, {40: 3.0}])
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    injection = MS_Injection({'SampleName': '2 blank runs'}, chromatogram, None, scans)

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(injection, out)

    _, sample_name = read_scans(out)
    assert sample_name.endswith('2_blank_runs')


def test_mzml_handles_a_missing_sample_name(tmp_path):
    '''The AIA fallback path synthesises None metadata, which must still produce a valid id.'''

    scans = make_scans([1000.0, 1500.0, 2000.0], [{40: 1.0}, {40: 2.0}, {40: 3.0}])
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    injection = MS_Injection({}, chromatogram, None, scans)

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(injection, out)

    round_tripped, sample_name = read_scans(out)
    pd.testing.assert_frame_equal(round_tripped, scans)
    assert sample_name == 'run'


def test_mzml_rejects_an_injection_without_scans(tmp_path):
    '''An injection carrying no scan matrix cannot be exported.'''

    chromatogram = np.array([[0.0, 0.1, 0.2, 0.3], [1.0, 1.0, 1.0, 1.0]])
    injection = MS_Injection({'SampleName': 'SMP-A1'}, chromatogram, None, None)

    with pytest.raises(ValueError, match='no scans'):
        write_injection_to_mzml(injection, tmp_path / 'out.mzML')


def test_write_sequence_to_mzml_writes_one_file_per_injection(ms_injection, tmp_path):
    '''Each injection in a sequence becomes its own round-trippable file.'''

    other = MS_Injection({'SampleName': 'SMP-A2'}, ms_injection.chromatogram, None,
                         ms_injection.scans)
    sequence = MS_Sequence({'sequence_name': 'SEQ'},
                           {'SMP-A1': ms_injection, 'SMP-A2': other})

    write_sequence_to_mzml(sequence, tmp_path)

    assert sorted(p.name for p in tmp_path.glob('*.mzML')) == ['SMP-A1.mzML', 'SMP-A2.mzML']
    scans, _ = read_scans(tmp_path / 'SMP-A2.mzML')
    pd.testing.assert_frame_equal(scans, ms_injection.scans)


# --- FID / ANDI-AIA netCDF ------------------------------------------------------------------

@pytest.fixture
def fid_injection():
    '''An FID_Injection over a synthetic, uniformly sampled chromatogram.'''

    return FID_Injection({'SampleName': 'FID-A1'}, make_fid_chromatogram(), 1.0)


def read_cdf(path):
    '''Reads a chromatogram back through pyGecko's own ANDI reader.

    The reader is a name-mangled private static; calling it directly is what makes this a
    round-trip test rather than a test of the writer against itself.
    '''

    return Agilent_FID_Parser._Agilent_FID_Parser__read_cdf_file(Path(path))


def test_cdf_round_trip_preserves_the_chromatogram(fid_injection, tmp_path):
    '''Writing then reading an FID injection reproduces the 2xN chromatogram.'''

    out = tmp_path / 'out.cdf'
    write_injection_to_cdf(fid_injection, out)

    round_tripped = read_cdf(out)

    assert round_tripped.shape == fid_injection.chromatogram.shape
    np.testing.assert_allclose(round_tripped, fid_injection.chromatogram, rtol=1e-9, atol=1e-9)


def test_cdf_writes_the_andi_variables(fid_injection, tmp_path):
    '''The file carries the ANDI/AIA variables and attributes, independent of pyGecko's reader.'''

    out = tmp_path / 'out.cdf'
    write_injection_to_cdf(fid_injection, out)

    dataset = nc.Dataset(out, 'r')
    try:
        assert 'ordinate_values' in dataset.variables
        interval = float(np.asarray(dataset.variables['actual_sampling_interval'][:]).item())
        delay = float(np.asarray(dataset.variables['actual_delay_time'][:]).item())
        intensities = np.asarray(dataset.variables['ordinate_values'][:], dtype=float)
        assert dataset.detector_unit == 'Arbitrary Intensity Units'
        assert dataset.retention_unit == 'Seconds'
        assert dataset.sample_name == 'FID-A1'
    finally:
        dataset.close()

    time = fid_injection.chromatogram[0]
    assert interval == pytest.approx(np.diff(time).mean() * 60.0)
    assert delay == pytest.approx(time[0] * 60.0)
    np.testing.assert_allclose(intensities, fid_injection.chromatogram[1])


def test_cdf_rejects_a_non_uniform_time_axis(tmp_path):
    '''ANDI reconstructs time from one interval, so a non-uniform axis cannot be represented.'''

    # Two different spacings spliced together: uniform enough to survive FID_Injection's
    # solvent-delay truncation, non-uniform where ANDI needs a single interval.
    time = np.concatenate([np.linspace(0.0, 5.0, 2000), np.linspace(5.01, 10.0, 500)])
    chromatogram = np.array([time, np.full(time.size, 5.0)])
    injection = FID_Injection({'SampleName': 'FID-A1'}, chromatogram, 1.0)

    with pytest.raises(ValueError, match='uniform'):
        write_injection_to_cdf(injection, tmp_path / 'out.cdf')


def test_write_sequence_to_cdf_writes_one_file_per_injection(fid_injection, tmp_path):
    '''Each FID injection in a sequence becomes its own round-trippable file.'''

    other = FID_Injection({'SampleName': 'FID-A2'}, make_fid_chromatogram(), 1.0)
    sequence = FID_Sequence({'sequence_name': 'SEQ'},
                            {'FID-A1': fid_injection, 'FID-A2': other})

    write_sequence_to_cdf(sequence, tmp_path)

    assert sorted(p.name for p in tmp_path.glob('*.cdf')) == ['FID-A1.cdf', 'FID-A2.cdf']
    np.testing.assert_allclose(read_cdf(tmp_path / 'FID-A2.cdf'), other.chromatogram,
                               rtol=1e-9, atol=1e-9)
