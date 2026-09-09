'''Round-trip tests for the export writers.

Every test here writes a file and reads it back with the reader pyGecko already ships, so the
success criterion is the readers' own output rather than a hand-checked byte layout. The export
is a nominal-mass reduction of the source, so round-trip fidelity is only ever claimed against
pyGecko's readers -- never against the original vendor file.
'''

from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd
import pytest

from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.injection.ms_injection import MS_Injection
from pygecko.gc_tools.sequence.fid_sequence import FID_Sequence
from pygecko.gc_tools.sequence.ms_sequence import MS_Sequence
from pygecko.parsers.agilent_fid_parser import Agilent_FID_Parser
from pygecko.parsers.file_readers import extract_scans_from_mzml
from pygecko.parsers.file_writers import (write_injection_to_cdf, write_injection_to_mzml,
                                          write_sequence_to_cdf, write_sequence_to_mzml)

from .conftest import fixture_path
from ..unit.conftest import make_fid_chromatogram, make_scans

pytestmark = pytest.mark.mzml

REAL_MZML = fixture_path('test_ri_calibration', 'FKB-FA-060-A1.mzML')


@pytest.fixture
def ms_injection():
    '''An MS_Injection over three scans, with the zero-filled gaps the readers produce.'''

    scans = make_scans(
        [89575.0, 89926.0, 90278.0],
        [{40: 1000.0, 41: 250.5}, {40: 900.0, 78: 12.0}, {41: 5.0, 78: 3000.25}],
    )
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    return MS_Injection({'SampleName': 'SMP-A1'}, chromatogram, None, scans)


def test_mzml_round_trip_preserves_scans(ms_injection, tmp_path):
    '''Writing then reading an injection reproduces the scans DataFrame exactly.'''

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    scans, sample_name = extract_scans_from_mzml(out)

    pd.testing.assert_frame_equal(scans, ms_injection.scans)
    assert sample_name == 'SMP-A1'


def test_mzml_round_trip_preserves_real_injection(tmp_path):
    '''A full 1878-scan injection read from a real mzML survives the write/read cycle.'''

    original, _ = extract_scans_from_mzml(Path(REAL_MZML))

    chromatogram = np.array([original.index / 60000, original.sum(axis=1)])
    injection = MS_Injection({'SampleName': 'FKB-FA-060-A1'}, chromatogram, None, original)

    out = tmp_path / 'FKB-FA-060-A1.mzML'
    write_injection_to_mzml(injection, out)
    round_tripped, sample_name = extract_scans_from_mzml(out)

    # check_dtype is off because the reader itself is dtype-inconsistent: pymzml returns float32
    # for two of this fixture's 241 m/z columns and float64 for the rest. The writer normalises
    # everything to float64, which is lossless, so the values still compare exactly.
    pd.testing.assert_frame_equal(round_tripped, original, check_dtype=False)
    assert sample_name == 'FKB-FA-060-A1'


def test_mzml_is_readable_by_an_independent_parser(ms_injection, tmp_path):
    '''pyteomics agrees with pymzml on the written file, so the output is not pymzml-specific.'''

    from pyteomics import mzml

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    with mzml.read(str(out)) as reader:
        spectra = list(reader)

    assert len(spectra) == len(ms_injection.scans)
    assert spectra[0]['ms level'] == 1
    # Zeros are dropped on write and restored by the reader's fillna, so the written arrays
    # only carry the m/z actually present in that scan.
    assert list(spectra[0]['m/z array']) == [40.0, 41.0]


def test_mzml_drops_zero_intensities(ms_injection, tmp_path):
    '''The dense zero-filled matrix is written sparsely.'''

    from pyteomics import mzml

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(ms_injection, out)

    with mzml.read(str(out)) as reader:
        lengths = [len(s['m/z array']) for s in reader]

    assert lengths == [2, 2, 2]
    assert ms_injection.scans.shape[1] == 3


def test_mzml_sanitises_a_sample_name_that_is_not_a_valid_xml_id(tmp_path):
    '''A sample name with spaces or a leading digit still produces a readable file.'''

    scans = make_scans([1000.0, 1500.0, 2000.0], [{40: 1.0}, {40: 2.0}, {40: 3.0}])
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    injection = MS_Injection({'SampleName': '2 blank runs'}, chromatogram, None, scans)

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(injection, out)

    _, sample_name = extract_scans_from_mzml(out)
    assert sample_name.endswith('2_blank_runs')


def test_mzml_handles_a_missing_sample_name(tmp_path):
    '''The AIA fallback path synthesises None metadata, which must still produce a valid id.'''

    scans = make_scans([1000.0, 1500.0, 2000.0], [{40: 1.0}, {40: 2.0}, {40: 3.0}])
    chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
    injection = MS_Injection({}, chromatogram, None, scans)

    out = tmp_path / 'out.mzML'
    write_injection_to_mzml(injection, out)

    round_tripped, sample_name = extract_scans_from_mzml(out)
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
    scans, _ = extract_scans_from_mzml(tmp_path / 'SMP-A2.mzML')
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
