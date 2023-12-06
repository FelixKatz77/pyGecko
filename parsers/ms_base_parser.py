import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from gc_tools.injection.ms_injection import MS_Injection
from parsers.msconvert_wraper import msconvert
from parsers.file_readers import extract_scans_from_mzml, extract_scans_from_mzxml


class MS_Base_Parser:

    @classmethod
    def load_injection(cls, raw_data_path: str, temp_dir: tempfile.TemporaryDirectory = None) -> MS_Injection:
        injection = MS_Base_Parser.initialize_injection(Path(raw_data_path), temp_dir=temp_dir)
        return injection

    @staticmethod
    def initialize_injection(path: Path,
                               temp_dir: tempfile.TemporaryDirectory = None) -> MS_Injection:
        '''
        Returns an MS_Injection object.

        Args:
            path (Path): Path to the raw directory.

        Returns:
            MS_Injection: An MS_Injection object
        '''

        scans = MS_Base_Parser.extract_scans_from_raw_data(path, temp_dir=temp_dir)
        chromatogram = MS_Base_Parser.construct_chromatogram_from_scans(scans)
        sample_name = path.stem
        injection = MS_Injection({'SampleName':sample_name}, chromatogram, None, scans)
        return injection


    @staticmethod
    def extract_scans_from_raw_data(raw_path: Path, temp_dir: tempfile.TemporaryDirectory = None) -> pd.DataFrame:

        '''
        Takes in the path to a raw file containing the scans of an injection, returns a DataFrame containing the scans.
        '''

        if not temp_dir:
            temp_dir = tempfile.TemporaryDirectory()
        if raw_path.suffix == '.mzML':
            try:
                scans_df = extract_scans_from_mzml(raw_path)
                return scans_df
            except KeyError as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise KeyError(f'Cannot extract scans from {raw_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)
        elif raw_path == '.mzXML':
            try:
                scans_df = extract_scans_from_mzxml(raw_path)
                return scans_df
            except KeyError as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise KeyError(f'Cannot extract scans from {raw_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)
        else:
            mzml_path = Path.joinpath(Path(temp_dir.name), raw_path.name).with_suffix('.mzML')
            try:
                msconvert([raw_path], temp_dir.name)
                scans_df = extract_scans_from_mzml(mzml_path)
                return scans_df
            except KeyError as error:
                raise KeyError(f'Cannot extract scans from {mzml_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)

    @staticmethod
    def construct_chromatogram_from_scans(scans: pd.DataFrame) -> np.ndarray:

        '''
        Returns a total ion chromatogram.

        Args:
            scans(pd.DataFrame): A DataFrame containing the scans of an injection.

        Returns:
            np.ndarray: The total ion chromatogram of the scans.
        '''

        chromatogram = np.array([scans.index, scans.sum(axis=1)])
        chromatogram[0] = chromatogram[0] / 60000
        return chromatogram