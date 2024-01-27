import numpy as np
from pathlib import Path
from pyteomics import mzxml
import pandas as pd
import pymzml
from pygecko.parsers.utilities import HiddenPrints


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



