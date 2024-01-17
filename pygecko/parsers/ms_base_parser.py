import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from pygecko.gc_tools import MS_Injection, RI_Calibration
from pygecko.parsers.msconvert_wraper import msconvert
from pygecko.parsers.file_readers import extract_scans_from_mzxml, extract_scans_from_mzml


class MS_Base_Parser:

    @classmethod
    def load_injection(cls, raw_data_path: str, temp_dir: tempfile.TemporaryDirectory = None) -> MS_Injection:

        '''
        Returns an MS_Injection object.

        Args:
            raw_data_path (str): Path to the raw data.
            temp_dir (tempfile.TemporaryDirectory): Temporary directory to store the converted mzML file. Defaults to None.

        Returns:
            MS_Injection: An MS_Injection object
        '''
        injection = MS_Base_Parser.initialize_injection(Path(raw_data_path), temp_dir=temp_dir)
        return injection

    @classmethod
    def load_ri_calibration(cls, raw_data_path: str, c_count: int, rt: float) -> RI_Calibration:

        '''
        Returns an RI_Calibration object.

        Args:
            raw_directory (str): Path to the directory containing the raw data.
            c_count (int): Number of carbon atoms for as specific alkane present in the standard.
            rt (float): Retention time of the alkane the c_count is provided for.

        Returns:
            RI_Calibration: An RI_Calibration object
        '''

        injection = MS_Base_Parser.initialize_injection(Path(raw_data_path))
        return RI_Calibration(injection, c_count, rt)

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