import numpy as np
from pathlib import Path
from pyteomics import mzxml
import pandas as pd
import pymzml
from pygecko.parsers.utilities import HiddenPrints
import netCDF4 as nc


def extract_scans_from_mzxml(mzxml_file: Path) -> pd.DataFrame:

    '''
    Takes in the path to a mzxml file containing the scans of an injection, returns a DataFrame containing the scans and
    the sample name.

    Args:
        mzxml_file (Path): Path to the mzxml file.

    Returns:
        tuple[pd.DataFrame, str]: DataFrame containing the scans and the sample name.
    '''

    scans = []
    with mzxml.read(str(mzxml_file)) as reader:
        sample_name = mzxml_file.stem
        for spectrum in reader:
            scan = {
                'retention_time': int(spectrum['retentionTime']*60000),
            }
            for m, i in zip(spectrum['m/z array'], spectrum['intensity array']):
                scan[round(m, 0)] = i
            scans.append(scan)

    df = pd.DataFrame(scans)
    df.fillna(0, inplace=True)
    df.set_index('retention_time', inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    return df, sample_name

def extract_scans_from_mzml(mzml_file: Path) -> (pd.DataFrame, str):

    '''
    Takes in the path to a mzml file containing the scans of an injection, returns a DataFrame containing the scans and
    the sample name.

    Args:
        mzml_file (Path): Path to the mzml file.

    Returns:
        tuple[pd.DataFrame, str]: DataFrame containing the scans and the sample name.
    '''
    with HiddenPrints():
        scans = []
        with pymzml.run.Reader(str(mzml_file)) as run:
            sample_name = run.info['run_id']
            for spectrum in run:
                retention_time = spectrum.scan_time[0] * 60000
                mzs = np.round(spectrum.mz).astype(int)
                intensities = spectrum.i

                scan = {'retention_time': retention_time}
                scan.update(dict(zip(mzs, intensities)))
                scans.append(scan)
                # scan = {
                #     'retention_time': spectrum.scan_time[0]*60000
                # }
                # for m, i in zip(spectrum.mz, spectrum.i):
                #     scan[round(m, 0)] = i
                # scans.append(scan)
        df = pd.DataFrame(scans).fillna(0).set_index('retention_time')
        df = df.reindex(sorted(df.columns), axis=1)
        return df, sample_name


def extract_scans_from_cdf(cdf_file: Path) -> tuple[pd.DataFrame, str]:
    '''
    Takes in the path to an AIA NetCDF (.cdf) file containing the scans of an injection,
    returns a DataFrame containing the scans and the sample name.

    Args:
        cdf_file (Path): Path to the .cdf file.

    Returns:
        tuple[pd.DataFrame, str]: DataFrame containing the scans and the sample name.
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

        num_scans = len(times)
        # One scatter over the flat ANDI arrays instead of a DataFrame per scan: the per-scan
        # groupby made this reader ~16x slower than the mzML path on a 1900-scan run. Each centroid
        # is binned to its nominal mass and the maximum intensity per (scan, mass) is kept.
        nominal_masses = np.round(np.asarray(masses)).astype(int)
        intensities = np.asarray(intensities, dtype=float)
        scan_lengths = np.diff(np.append(np.asarray(scan_index), len(nominal_masses)))
        rows = np.repeat(np.arange(num_scans), scan_lengths)
        columns, column_indices = np.unique(nominal_masses, return_inverse=True)
        matrix = np.zeros((num_scans, len(columns)))
        np.maximum.at(matrix, (rows, column_indices), intensities)

        # AIA NetCDF time is usually in seconds. Multiply by 1000 for milliseconds.
        retention_times = np.asarray(times, dtype=float) * 1000

    finally:
        dataset.close()  # Ensure the file is closed even if an error occurs

    # Format the DataFrame exactly like the mzML output
    df = pd.DataFrame(matrix, columns=columns, index=pd.Index(retention_times, name='retention_time'))

    return df, sample_name