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
        scans = []

        for i in range(num_scans):
            start_idx = scan_index[i]
            # The end index is either the start of the next scan, or the end of the array
            end_idx = scan_index[i + 1] if i < num_scans - 1 else len(masses)

            # Extract this specific scan's data
            scan_masses = np.round(masses[start_idx:end_idx]).astype(int)
            scan_intensities = intensities[start_idx:end_idx]

            # AIA NetCDF time is usually in seconds. Multiply by 1000 for milliseconds.
            retention_time = times[i] * 1000

            # scan = {'retention_time': retention_time}
            # scan.update(dict(zip(scan_masses, scan_intensities)))
            # scans.append(scan)

            scan = {'retention_time': retention_time}

            # Create a temporary DataFrame for the scan to handle duplicate rounded masses
            scan_df = pd.DataFrame({'mz': scan_masses, 'intensity': scan_intensities})

            # Group by the nominal mass and take the maximum intensity in that bin
            binned_scan = scan_df.groupby('mz')['intensity'].max().to_dict()

            scan.update(binned_scan)
            scans.append(scan)

    finally:
        dataset.close()  # Ensure the file is closed even if an error occurs

    # Format the DataFrame exactly like the mzML output
    df = pd.DataFrame(scans).fillna(0).set_index('retention_time')
    df = df.reindex(sorted(df.columns), axis=1)

    return df, sample_name