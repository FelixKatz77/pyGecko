import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from pygecko.gc_tools import MS_Injection, RI_Calibration, MS_Sequence
from pygecko.parsers.msconvert_wraper import msconvert
from pygecko.parsers.file_readers import extract_scans_from_mzxml, extract_scans_from_mzml


class MS_Base_Parser:

    @staticmethod
    def load_sequence(raw_directory: Path | str, pos: bool = False) -> MS_Sequence:

        '''
        Returns an MS_Sequence object.

        Args:
            raw_directory (Path|str): Path to a directory containing raw data.

        Returns:
            MS_Sequence: An MS_Sequence object.
        '''

        print('Loading GC-MS sequence...')
        raw_directory = Path(raw_directory)
        supported_formats = ['.D', '.mzML', '.mzXML']
        raw_files = []
        for file_format in supported_formats:
            raw_files.extend(raw_directory.glob(f'*{file_format}'))
        injections = {}
        for raw_file in raw_files:
            injection = MS_Base_Parser.load_injection(raw_file, pos=pos)
            injections[injection.sample_name] = injection
        print(f'Sequence loaded with {len(injections)} injections.')
        return MS_Sequence({}, injections)

    @staticmethod
    def load_injection(raw_data_path: str|Path, pos:bool=False, temp_dir: tempfile.TemporaryDirectory = None) -> MS_Injection:

        '''
        Returns an MS_Injection object.

        Args:
            raw_data_path (str): Path to the raw data.
            temp_dir (tempfile.TemporaryDirectory): Temporary directory to store the converted mzML file. Defaults to None.

        Returns:
            MS_Injection: An MS_Injection object
        '''
        injection = MS_Base_Parser.initialize_injection(Path(raw_data_path), temp_dir=temp_dir, pos=pos)
        return injection

    @staticmethod
    def load_ri_calibration(raw_data_path: str, c_count: int, rt: float) -> RI_Calibration:

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
    def initialize_injection(path: Path, pos:bool=False,
                               temp_dir: tempfile.TemporaryDirectory = None) -> MS_Injection:
        '''
        Returns an MS_Injection object.

        Args:
            path (Path): Path to the raw directory.

        Returns:
            MS_Injection: An MS_Injection object
        '''

        scans, sample_name = MS_Base_Parser.extract_scans_from_raw_data(path, temp_dir=temp_dir)
        chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
        injection = MS_Injection({'SampleName':sample_name}, chromatogram, None, scans, pos=pos)
        return injection


    @staticmethod
    def extract_scans_from_raw_data(raw_path: Path, temp_dir: tempfile.TemporaryDirectory = None) -> (pd.DataFrame, str):

        '''
        Takes in the path to a raw file containing the scans of an injection, returns a DataFrame containing the scans
        and the injection's sample name.
        '''


        if raw_path.suffix == '.mzML':
            try:
                scans_df, sample_name = extract_scans_from_mzml(raw_path)
                return scans_df, sample_name
            except KeyError as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise KeyError(f'Cannot extract scans from {raw_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)
        elif raw_path.suffix == '.mzXML':
            try:
                scans_df, sample_name = extract_scans_from_mzxml(raw_path)
                return scans_df, sample_name
            except KeyError as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise KeyError(f'Cannot extract scans from {raw_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)
        else:
            if not temp_dir:
                temp_dir = tempfile.TemporaryDirectory()
            mzml_path = Path.joinpath(Path(temp_dir.name), raw_path.name).with_suffix('.mzML')
            try:
                msconvert([raw_path], temp_dir.name)
                scans_df, sample_name = extract_scans_from_mzml(mzml_path)
                return scans_df, sample_name
            except KeyError as error:
                raise KeyError(f'Cannot extract scans from {mzml_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)

