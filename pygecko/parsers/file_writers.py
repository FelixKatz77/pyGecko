import re
from pathlib import Path

import netCDF4 as nc
import numpy as np
from psims.mzml import MzMLWriter

from pygecko import __version__
from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.injection.ms_injection import MS_Injection
from pygecko.gc_tools.sequence.fid_sequence import FID_Sequence
from pygecko.gc_tools.sequence.ms_sequence import MS_Sequence

def _xml_id(name: str) -> str:

    '''
    Returns a sample name reduced to a valid XML ID.

    The mzML schema types the run's id attribute as xsd:ID, so it may not contain whitespace and
    may not begin with a digit. Sample names out of the vendor metadata are unconstrained, so they
    are sanitised rather than trusted.

    Args:
        name (str): Sample name, possibly None.

    Returns:
        str: A valid XML ID derived from the name.
    '''

    if not name:
        return 'run'
    sanitised = re.sub(r'[^\w.-]', '_', name, flags=re.ASCII)
    if not re.match(r'[A-Za-z_]', sanitised):
        sanitised = f'run_{sanitised}'
    return sanitised


def write_injection_to_mzml(injection: MS_Injection, path: Path|str) -> None:

    '''
    Writes an MS injection to an indexed mzML file.

    The spectra are the centroids exactly as the reader delivered them: unrounded m/z and the
    intensities in their source dtype, so the file holds what the vendor file held. Polarity,
    acquisition start time and instrument are written only when the source provided them; a term
    the source lacked is left out rather than defaulted. An injection that carries no raw scans
    (one saved before they were kept) is written from its nominal-mass matrix instead, with the
    binning declared in the file's processing history, and the zeros that square off the matrix
    are dropped on write.

    Args:
        injection (MS_Injection): Injection to write.
        path (Path|str): Path of the mzML file to write.

    Raises:
        ValueError: If the injection holds no scan matrix.
    '''

    if injection.scans is None:
        raise ValueError(f'Cannot write {injection.sample_name} to mzML: the injection holds no '
                         f'scans.')

    if injection.raw_scans is not None:
        spectra = injection.raw_scans.spectra()
        processing = ['Conversion to mzML']
    else:
        spectra = _matrix_spectra(injection.scans)
        processing = ['Conversion to mzML', {'nominal mass binning': None}]
    run_id = _xml_id(injection.sample_name)
    start_time = injection.acq_time.isoformat() if injection.acq_time else None
    instrument = [injection.instrument_name] if injection.instrument_name else []

    with MzMLWriter(str(path)) as writer:
        writer.controlled_vocabularies()
        writer.file_description(['MS1 spectrum', 'centroid spectrum'])
        writer.software_list([{'id': 'pygecko', 'version': __version__,
                               'params': ['custom unreleased software tool']}])
        writer.instrument_configuration_list(
            [writer.InstrumentConfiguration(id='IC1', component_list=[], params=instrument)])
        writer.data_processing_list([writer.DataProcessing(
            [writer.ProcessingMethod(order=1, software_reference='pygecko', params=processing)],
            id='DP1')])

        with writer.run(id=run_id, instrument_configuration='IC1', start_time=start_time):
            with writer.spectrum_list(count=len(injection.scans)):
                for index, (rt_ms, mzs, intensities) in enumerate(spectra, start=1):
                    writer.write_spectrum(
                        mzs, intensities,
                        id=f'scan={index}',
                        polarity=injection.polarity,
                        centroided=True,
                        scan_start_time=rt_ms / 60000,
                        encoding={'m/z array': np.float64,
                                  'intensity array': intensities.dtype.type},
                        params=[{'ms level': 1}, {'total ion current': intensities.sum()}])
            with writer.chromatogram_list(count=1):
                writer.write_chromatogram(
                    injection.chromatogram[0], injection.chromatogram[1],
                    id='TIC', chromatogram_type='total ion current chromatogram',
                    encoding=np.float64)


def _matrix_spectra(scans):

    '''
    Yields the retention time in milliseconds, m/z array and intensity array of every row of a
    nominal-mass matrix, without the zeros that square it off.
    '''

    mzs = scans.columns.to_numpy(dtype=np.float64)
    for rt_ms, row in zip(scans.index, scans.to_numpy(dtype=np.float64)):
        present = row != 0
        yield rt_ms, mzs[present], row[present]


def write_sequence_to_mzml(sequence: MS_Sequence, directory: Path|str) -> None:

    '''
    Writes every injection of an MS sequence to its own mzML file, named after the sample.

    Args:
        sequence (MS_Sequence): Sequence to write.
        directory (Path|str): Directory to write the files into. Created if it does not exist.

    '''

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for injection in sequence.injections.values():
        write_injection_to_mzml(injection, directory / f'{injection.sample_name}.mzML')


def write_injection_to_cdf(injection: FID_Injection, path: Path|str) -> None:

    '''
    Writes an FID injection to an ANDI/AIA netCDF file.

    FID data is written as netCDF rather than mzML deliberately. The PSI-MS controlled vocabulary
    holds no term for a flame ionization detector and none of its chromatogram types describes one,
    so an mzML export would be schema-valid but semantically wrong. ANDI/AIA (ASTM E1947/E1948) is
    the chromatography standard for a detector trace, and it is a format pyGecko already reads.

    The format reconstructs the time axis from a single sampling interval, so the chromatogram's
    time axis must be uniformly sampled.

    Args:
        injection (FID_Injection): Injection to write.
        path (Path|str): Path of the .cdf file to write.

    Raises:
        ValueError: If the chromatogram's time axis is not uniformly sampled.
    '''

    time_min, intensities = injection.chromatogram
    intervals = np.diff(time_min)
    if not np.allclose(intervals, intervals[0], rtol=1e-6, atol=1e-12):
        raise ValueError(f'Cannot write {injection.sample_name} to netCDF: ANDI/AIA reconstructs '
                         f'the time axis from a single sampling interval, so the chromatogram must '
                         f'be uniform.')

    dataset = nc.Dataset(path, 'w', format='NETCDF3_CLASSIC')
    try:
        dataset.dataset_completeness = 'C1'
        dataset.aia_template_revision = '1.0'
        dataset.netcdf_revision = '2.3.2'
        dataset.sample_name = injection.sample_name or ''
        dataset.detector_unit = 'Arbitrary Intensity Units'
        dataset.retention_unit = 'Seconds'
        dataset.detector_name = 'FID'

        dataset.createDimension('point_number', len(intensities))

        ordinates = dataset.createVariable('ordinate_values', 'f8', ('point_number',))
        ordinates[:] = intensities

        interval = dataset.createVariable('actual_sampling_interval', 'f8')
        interval.assignValue(float(intervals.mean() * 60.0))

        delay = dataset.createVariable('actual_delay_time', 'f8')
        delay.assignValue(float(time_min[0] * 60.0))

        run_time = dataset.createVariable('actual_run_time_length', 'f8')
        run_time.assignValue(float((time_min[-1] - time_min[0]) * 60.0))
    finally:
        dataset.close()


def write_sequence_to_cdf(sequence: FID_Sequence, directory: Path|str) -> None:

    '''
    Writes every injection of an FID sequence to its own ANDI/AIA netCDF file, named after the
    sample.

    Args:
        sequence (FID_Sequence): Sequence to write.
        directory (Path|str): Directory to write the files into. Created if it does not exist.
    '''

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for injection in sequence.injections.values():
        write_injection_to_cdf(injection, directory / f'{injection.sample_name}.cdf')
