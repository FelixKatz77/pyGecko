import tempfile
from pathlib import Path

import numpy as np

from pygecko.gc_tools import MS_Injection, RI_Calibration, MS_Sequence
from pygecko.gc_tools.injection.raw_scans import Raw_Scans
from pygecko.parsers.msconvert_wraper import msconvert
from pygecko.parsers.file_readers import extract_scans_from_mzxml, extract_scans_from_mzml, extract_scans_from_cdf
from typing import Iterable, Optional


class MS_Base_Parser:

    @staticmethod
    def load_sequence(
            raw_directory: Path | str,
            pos: bool = False,
            sample_filter: Optional[Iterable[str]] = None,
    ) -> MS_Sequence:
        """Returns an MS_Sequence object.

        Args:
            raw_directory: Path to a directory containing raw data.
            pos: Indicates if plate position is given in the injection names.
            sample_filter: Optional iterable of allowed sample names. Injections
                whose ``sample_name`` is not in this set are silently dropped.
                Useful to exclude cleaning/conditioning runs that may otherwise
                be ingested as real injections. Defaults to None (keep all).

        Returns:
            An MS_Sequence object.
        """

        print('Loading GC-MS sequence...')
        raw_directory = Path(raw_directory)
        supported_formats = ['.D', '.mzML', '.mzXML', 'cdf', '.CDF']
        raw_files = []
        for file_format in supported_formats:
            raw_files.extend(raw_directory.glob(f'*{file_format}'))

        allowed_names = set(sample_filter) if sample_filter is not None else None

        injections = {}
        for raw_file in raw_files:
            if raw_file.suffix.lower() == '.cdf' and '_spectra' not in raw_file.name:
                continue  # Skip FID or TIC-only CDF files
            injection = MS_Base_Parser.load_injection(raw_file, pos=pos)
            if allowed_names is not None and injection.sample_name not in allowed_names:
                continue
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

        raw_scans, metadata = MS_Base_Parser.extract_scans_from_raw_data(path, temp_dir=temp_dir)
        scans = raw_scans.to_nominal_matrix()
        chromatogram = np.array([scans.index / 60000, scans.sum(axis=1)])
        injection = MS_Injection(metadata, chromatogram, None, scans, pos=pos, raw_scans=raw_scans)
        injection.record_step('MS_Base_Parser.load_injection',
                              {'raw_data_path': str(path), 'pos': pos})
        return injection


    @staticmethod
    def extract_scans_from_raw_data(raw_path: Path, temp_dir: tempfile.TemporaryDirectory = None) -> tuple[Raw_Scans, dict]:

        '''
        Takes in the path to a raw file containing the scans of an injection, returns the centroids as read and the
        injection's metadata (SampleName, AcqTime, Polarity, InstrumentName).
        '''


        if raw_path.suffix == '.mzML':
            try:
                raw_scans, metadata = extract_scans_from_mzml(raw_path)
                return raw_scans, metadata
            except KeyError as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise KeyError(f'Cannot extract scans from {raw_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)
        elif raw_path.suffix == '.mzXML':
            try:
                raw_scans, metadata = extract_scans_from_mzxml(raw_path)
                return raw_scans, metadata
            except KeyError as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise KeyError(f'Cannot extract scans from {raw_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)
        elif raw_path.suffix.lower() == '.cdf':
            try:
                raw_scans, metadata = extract_scans_from_cdf(raw_path)
                return raw_scans, metadata
            except Exception as error:
                print(f'Cannot extract scans from {raw_path.name}: {error}')
                raise RuntimeError(f'Cannot extract scans from {raw_path.name}: {error}')
        else:
            if not temp_dir:
                temp_dir = tempfile.TemporaryDirectory()
            mzml_path = Path.joinpath(Path(temp_dir.name), raw_path.name).with_suffix('.mzML')
            try:
                msconvert([raw_path], temp_dir.name)
                raw_scans, metadata = extract_scans_from_mzml(mzml_path)
                return raw_scans, metadata
            except KeyError as error:
                raise KeyError(f'Cannot extract scans from {mzml_path.name}: {error}')
            except FileNotFoundError as error:
                raise FileNotFoundError(error)

