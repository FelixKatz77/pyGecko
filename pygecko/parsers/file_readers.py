from datetime import datetime
import numpy as np
from pathlib import Path
from pyteomics import mzxml
import pymzml
from pygecko.gc_tools.injection.raw_scans import Raw_Scans
from pygecko.parsers.utilities import HiddenPrints
import netCDF4 as nc


def _flatten(retention_times: list, mzs: list[np.ndarray], intensities: list[np.ndarray]) -> Raw_Scans:

    '''
    Takes in one retention time, m/z array and intensity array per scan, returns them as Raw_Scans.
    '''

    scan_index = np.cumsum([0] + [len(mz) for mz in mzs[:-1]])
    return Raw_Scans(retention_times, scan_index, np.concatenate(mzs), np.concatenate(intensities))


def _parse_iso_timestamp(timestamp: str|None) -> datetime|None:

    '''
    Takes in an ISO 8601 timestamp as mzML writes it, returns a datetime. Python 3.10's
    fromisoformat does not accept the Z suffix msConvert emits for UTC.
    '''

    if not timestamp:
        return None
    return datetime.fromisoformat(timestamp.replace('Z', '+00:00'))


def extract_scans_from_mzxml(mzxml_file: Path) -> tuple[Raw_Scans, dict]:

    '''
    Takes in the path to a mzxml file containing the scans of an injection, returns the centroids as
    read and the injection's metadata.

    mzXML carries no run start time or instrument description that pyteomics exposes, so only the
    sample name (the file stem) and the polarity are filled in.

    Args:
        mzxml_file (Path): Path to the mzxml file.

    Returns:
        tuple[Raw_Scans, dict]: The centroids and the metadata (SampleName, AcqTime, Polarity,
        InstrumentName).
    '''

    retention_times, mzs, intensities = [], [], []
    polarity = None
    with mzxml.read(str(mzxml_file)) as reader:
        for spectrum in reader:
            if polarity is None:
                polarity = {'+': 'positive', '-': 'negative'}.get(spectrum.get('polarity'))
            retention_times.append(spectrum['retentionTime'] * 60000)
            mzs.append(spectrum['m/z array'])
            intensities.append(spectrum['intensity array'])
    metadata = {'SampleName': mzxml_file.stem, 'AcqTime': None, 'Polarity': polarity,
                'InstrumentName': None}
    return _flatten(retention_times, mzs, intensities), metadata


def extract_scans_from_mzml(mzml_file: Path) -> tuple[Raw_Scans, dict]:

    '''
    Takes in the path to a mzml file containing the scans of an injection, returns the centroids as
    read and the injection's metadata.

    The polarity is taken from the first spectrum: a GC-MS run is acquired in one polarity. The
    instrument is the first parameter of the instrument configuration: the instrument-model
    cvParam when the file has one (msConvert writes it), or the userParam pyGecko's own writer
    emits for a free-text instrument name.

    Args:
        mzml_file (Path): Path to the mzml file.

    Returns:
        tuple[Raw_Scans, dict]: The centroids and the metadata (SampleName, AcqTime, Polarity,
        InstrumentName).
    '''

    retention_times, mzs, intensities = [], [], []
    polarity = None
    with HiddenPrints():
        with pymzml.run.Reader(str(mzml_file)) as run:
            for spectrum in run:
                if polarity is None:
                    if 'MS:1000130' in spectrum:
                        polarity = 'positive'
                    elif 'MS:1000129' in spectrum:
                        polarity = 'negative'
                retention_times.append(spectrum.scan_time[0] * 60000)
                mzs.append(spectrum.mz)
                intensities.append(spectrum.i)
            configuration = run.info['instrument_configuration_list_element'].find(
                '{*}instrumentConfiguration')
            param = configuration.find('{*}cvParam') if configuration is not None else None
            if param is None and configuration is not None:
                param = configuration.find('{*}userParam')
            instrument_name = param.attrib['name'] if param is not None else None
            metadata = {'SampleName': run.info['run_id'],
                        'AcqTime': _parse_iso_timestamp(run.info.get('start_time')),
                        'Polarity': polarity, 'InstrumentName': instrument_name}
    return _flatten(retention_times, mzs, intensities), metadata


def extract_scans_from_cdf(cdf_file: Path) -> tuple[Raw_Scans, dict]:
    '''
    Takes in the path to an AIA NetCDF (.cdf) file containing the scans of an injection,
    returns the centroids as read and the injection's metadata.

    Args:
        cdf_file (Path): Path to the .cdf file.

    Returns:
        tuple[Raw_Scans, dict]: The centroids and the metadata (SampleName, AcqTime, Polarity,
        InstrumentName).
    '''
    # Extract sample name (Assuming OpenLab format: SampleName_MS1Front...)
    sample_name = cdf_file.name.split('_')[0]

    dataset = nc.Dataset(cdf_file, 'r')

    try:
        # Read the standard ANDI-MS arrays
        times = dataset.variables['scan_acquisition_time'][:]  # Typically in seconds
        masses = dataset.variables['mass_values'][:]
        intensities = dataset.variables['intensity_values'][:]
        scan_index = dataset.variables['scan_index'][:]

        # AIA NetCDF time is usually in seconds. Multiply by 1000 for milliseconds.
        raw = Raw_Scans(np.asarray(times, dtype=float) * 1000, np.asarray(scan_index),
                        np.asarray(masses), np.asarray(intensities))

        # ANDI-MS global attributes; both are optional in the standard.
        polarity = getattr(dataset, 'test_ionization_polarity', '').lower()
        polarity = polarity.split()[0] if polarity.startswith(('positive', 'negative')) else None
        timestamp = getattr(dataset, 'experiment_date_time_stamp', None)  # YYYYMMDDhhmmss+hhmm
        acq_time = datetime.strptime(timestamp, '%Y%m%d%H%M%S%z') if timestamp else None

    finally:
        dataset.close()  # Ensure the file is closed even if an error occurs

    metadata = {'SampleName': sample_name, 'AcqTime': acq_time, 'Polarity': polarity,
                'InstrumentName': None}
    return raw, metadata
