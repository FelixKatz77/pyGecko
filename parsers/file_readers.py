from pathlib import Path
from pyteomics import mzxml, mzml
import pandas as pd

def extract_scans_from_mzxml(mzxml_file: Path) -> pd.DataFrame:

    '''
    Takes in the path to a mzxml file containing the scans of an injection, returns a DataFrame containing the scans.

    Args:
        mzxml_file (Path): Path to the mzxml file.

    Returns:
        pd.DataFrame: DataFrame containing the scans.
    '''

    scans = []
    with mzxml.read(str(mzxml_file)) as reader:
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
    return df

def extract_scans_from_mzml(mzml_file: Path) -> pd.DataFrame:

    '''
    Takes in the path to a mzml file containing the scans of an injection, returns a DataFrame containing the scans.

    Args:
        mzml_file (Path): Path to the mzml file.

    Returns:
        pd.DataFrame: DataFrame containing the scans.
    '''

    scans = []
    with mzml.read(str(mzml_file)) as reader:
        for spectrum in reader:
            scan = {
                'retention_time': int(spectrum['scanList']['scan'][0]['scan start time']*60000),
            }
            for m, i in zip(spectrum['m/z array'], spectrum['intensity array']):
                scan[round(m, 0)] = i
            scans.append(scan)

    df = pd.DataFrame(scans)
    df.fillna(0, inplace=True)
    df.set_index('retention_time', inplace=True)
    df = df.reindex(sorted(df.columns), axis=1)
    return df

if __name__ == "__main__":
    mzml_file = "C:/Users/felix/Downloads/25.11.2021_Bame Std_3.mzML"  # Replace with the path to your mzXML file
    data_df = extract_scans_from_mzml(mzml_file)
    print(data_df)