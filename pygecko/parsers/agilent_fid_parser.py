import xml.etree.ElementTree as ET
import re
import numpy as np
from pathlib import Path
from datetime import datetime

from pygecko.parsers.fid_base_parser import FID_Base_Parser
from pygecko.gc_tools import FID_Sequence, FID_Injection, RI_Calibration
import netCDF4 as nc
from typing import Iterable, Optional



class Agilent_FID_Parser:

    '''
    A class wrapping functions to parse raw data from an Agilent GC/FID into FID_Sequence objects.
    '''

    @classmethod
    def load_sequence(cls, raw_directory:str, solvent_delay:float|int, pos:bool=False, file_source: str='csv', sample_filter: Optional[Iterable[str]] = None) -> FID_Sequence:
        '''
        Returns an FID_Sequence object.

        Args:
            raw_directory (str): Path to the directory containing the corresponding xy-files and the acaml-file.
            solvent_delay (float): Retention time of the solvent peak in minutes.
            pos (bool, optional): Indicates if plate position is given in the injection names. Defaults to False.
            file_source (str, optional): Locating FID signal files. Defaults to 'csv' for legacy .xy/CSV in root, 'cdf' for column split system.
            Optional iterable of allowed ``SampleName`` values. Injections whose ``SampleName`` is not in this set are silently skipped. Useful for excluding cleaning/conditioning runs that
                are tagged ``SampleType='Sample'`` rather than ``'Blank'``. Defaults to None (no extra filtering).

        Returns:
            FID_Sequence: An FID_Sequence object
        '''

        print('Loading GC-FID sequence...')
        acaml_files = list(Path(raw_directory).glob('*.acaml'))
        if acaml_files:
            sequence_metadata, injections = cls.__load_sequence_data(
                acaml_files[0],
                raw_directory,
                solvent_delay,
                pos=pos,
                file_source=file_source,
                sample_filter=sample_filter,
            )
        else:
            # No acaml (sequence-level OpenLab metadata) in the result folder.
            # This happens with incomplete/partial exports. Fall back to
            # enumerating the FID injections straight from the AIA/*_FID1A.cdf
            # files, deriving the sample name from each filename.
            print(
                f'No .acaml file found in {raw_directory}; falling back to '
                f'enumerating FID injections from AIA/*_FID1A.cdf.'
            )
            sequence_metadata, injections = cls.__load_sequence_data_from_aia(
                raw_directory,
                solvent_delay,
                pos=pos,
                file_source=file_source,
                sample_filter=sample_filter,
            )
        print(f'Sequence loaded with {len(injections)} injections.')
        return FID_Sequence(sequence_metadata, injections)

    @classmethod
    def load_ri_calibration(cls, raw_directory:str, solvent_delay, c_count:int, rt:float, file_source: str = 'csv',) -> RI_Calibration:

        '''
        Returns an RI_Calibration object.

        Args:
            raw_directory (str): Path to the directory containing the corresponding xy-file and the acaml-file.
            solvent_delay (float): Retention time of the solvent peak in minutes.
            c_count (int): Carbon count of the alkane the retention time is provided for.
            rt (float): Retention time of the alkane the c_count is provided for.
            file_source: Strategy for locating the FID signal file. See ``load_sequence``. Defaults to ``'csv'``.

        Returns:
            RI_Calibration: An RI_Calibration object.
        '''

        acaml_file = list(Path(raw_directory).glob('*.acaml'))[0]
        injection = cls.__load_injection_data(acaml_file, raw_directory, solvent_delay, file_source=file_source)
        return RI_Calibration(injection, c_count, rt)

    @classmethod
    def load_injection(cls, raw_directory:str, solvent_delay:float|int, file_source: str = 'csv',) -> FID_Injection:

        '''
        Returns an FID_Injection object.

        Args:
            raw_directory (str): Path to the directory containing the corresponding xy-file and the acaml-file.
            solvent_delay (float): Retention time of the solvent peak in minutes.
            file_source: Strategy for locating the FID signal file. See ``load_sequence``. Defaults to ``'csv'``.

        Returns:
            FID_Injection: An FID_Injection object.
        '''

        acaml_file = list(Path(raw_directory).glob('*.acaml'))[0]
        injection = cls.__load_injection_data(acaml_file, raw_directory, solvent_delay, file_source=file_source)
        return injection

    @staticmethod
    def __load_injection_data(acaml_file:Path, raw_directory:str, solvent_delay:float|int, file_source: str = 'csv') -> FID_Injection:

        '''
        Returns an FID_Injection object.

        Args:
            acaml_file (Path): Path to an acaml-file from an Agilent GC.
            raw_directory (str): Path to the directory containing the corresponding xy-file and the acaml-file.
            solvent_delay (float): Retention time  of the solvent peak in minutes.
            file_source: Strategy for locating the FID signal file. See ``load_sequence``. Defaults to ``'csv'``.


        Returns:
            FID_Injection: An FID_Injection object.
        '''

        root = Agilent_FID_Parser.__get_root(acaml_file)
        metadata_array = list(root.iter(f'ArrayOfInjectionMetaData'))[0]
        metadata_element = list(metadata_array.iter('InjectionMetaData'))[0]
        injection_metadata = Agilent_FID_Parser.__get_injection_metadata(metadata_element)
        resolved_source = Agilent_FID_Parser.__resolve_file_source(raw_directory, file_source=file_source)
        xy_array = Agilent_FID_Parser.__get_fid_data_for_injection(
            raw_directory, injection_metadata, resolved_source,
        )
        if xy_array is None and resolved_source == 'csv':
            # Legacy fallback: directory may contain a single signal file with
            # a name unrelated to the SampleName attribute. Mirrors the original
            # __load_injection_data behavior. The two-extension list relies on
            # filesystem case-insensitivity (Windows) — listing more variants
            # would double-count the same file.
            raw_path = Path(raw_directory)
            candidates = list(raw_path.glob('*.xy')) + list(raw_path.glob('*.CSV'))
            if len(candidates) == 1:
                xy_array = FID_Base_Parser.read_xy_array(candidates[0])

        if xy_array is None:
            raise ValueError(
                f'Could not locate FID data for injection '
                f'{injection_metadata.get("SampleName")} in {raw_directory} '
                f'(file_source={resolved_source}).'
            )
        injection = FID_Injection(injection_metadata, xy_array, solvent_delay)
        return injection


    @staticmethod
    def __load_sequence_data(acaml_file:Path, xy_directory:str, solvent_delay:float|int, pos:bool=False, file_source: str = 'csv',
        sample_filter: Optional[Iterable[str]] = None,) -> (dict, dict[str:FID_Injection]):

        '''
        Returns a dict containing the sequence metadata and a dict containing the injections.

        Args:
            acaml_file (Path): Path to an acaml-file from an Agilent GC
            xy_directory (str): Path to the directory containing the corresponding xy-files.
            solvent_delay (float): Retention time  of the solvent peak in minutes.
            pos (bool, optional): Indicates if plate position is given in the injection names. Defaults to False.
            file_source: Strategy for locating per-injection FID signal files.
            sample_filter: Optional iterable of allowed ``SampleName`` values

        Returns:
            tuple[dict, dict[str:FID_Injection]]: A tuple containing a dict containing the sequence metadata and
            a dict containing the injections.
        '''

        root = Agilent_FID_Parser.__get_root(acaml_file)
        xlmns = Agilent_FID_Parser.__get_xmlns(root)
        sequence_metadata = Agilent_FID_Parser.__get_sequence_metadata(root, xlmns)
        injections_metadata = Agilent_FID_Parser.__get_injections_metadata(root, sample_filter=sample_filter)
        xy_arrays = Agilent_FID_Parser.__get_xy_arrays(
            xy_directory, injections_metadata, file_source=file_source)
        injections = Agilent_FID_Parser.__initialize_injections(injections_metadata, xy_arrays, solvent_delay, pos=pos)
        return sequence_metadata, injections

    @staticmethod
    def __load_sequence_data_from_aia(raw_directory: str, solvent_delay: float | int, pos: bool = False,
                                      file_source: str = 'cdf',
                                      sample_filter: Optional[Iterable[str]] = None) -> (dict, dict[str:FID_Injection]):

        '''
        Returns sequence metadata and injections enumerated from the AIA/*_FID1A.cdf files.

        Fallback for result folders that lack an acaml (the sequence-level OpenLab metadata file), e.g. from an
        incomplete or partial export. Instead of reading the injection list from the acaml, each FID injection is
        discovered directly from an ``AIA/<sample>_FID1A.cdf`` file, with the sample name taken from the filename
        (the part preceding ``_FID1A``). The retention-time axis is reconstructed from the cdf exactly as in the
        acaml-driven path. No acaml means no acquisition-method/instrument/timestamp metadata is available, so those
        fields are left as None.

        Args:
            raw_directory (str): Path to the result directory; must contain an ``AIA`` subdirectory with the
                ``*_FID1A.cdf`` traces.
            solvent_delay (float): Retention time of the solvent peak in minutes.
            pos (bool, optional): Indicates if plate position is given in the injection names. Defaults to False.
            file_source (str): FID ingestion strategy. Only ``'cdf'`` (or ``'auto'`` resolving to ``'cdf'``) is
                supported here, since enumeration relies on the AIA cdf files. Defaults to ``'cdf'``.
            sample_filter: Optional iterable of allowed sample names. Injections whose name is not in this set are
                skipped. Defaults to None (no extra filtering).

        Returns:
            tuple[dict, dict[str:FID_Injection]]: A tuple containing a dict with the sequence metadata and a dict
            containing the injections keyed by sample name.

        Raises:
            FileNotFoundError: If the ``AIA`` subdirectory or any ``*_FID1A.cdf`` files are missing, or if
                ``file_source`` does not resolve to ``'cdf'`` (the only mode enumeration supports without an acaml).
        '''

        resolved_source = Agilent_FID_Parser.__resolve_file_source(raw_directory, file_source)
        if resolved_source != 'cdf':
            raise FileNotFoundError(
                f'No .acaml file found in {raw_directory} and file_source={file_source!r} does not resolve to '
                f"'cdf'. Without an acaml the FID injections can only be enumerated from AIA/*_FID1A.cdf files; "
                f"pass file_source='cdf' (or restore the .acaml)."
            )
        aia_path = Path(raw_directory) / 'AIA'
        if not aia_path.exists():
            raise FileNotFoundError(
                f'No .acaml file and no AIA subdirectory found in {raw_directory}; cannot enumerate FID injections.'
            )
        cdf_files = sorted(aia_path.glob('*_FID1A.cdf'))
        if not cdf_files:
            raise FileNotFoundError(
                f'No .acaml file and no *_FID1A.cdf files found in {aia_path}; cannot enumerate FID injections.'
            )

        allowed_names = set(sample_filter) if sample_filter is not None else None
        injections = {}
        for cdf_path in cdf_files:
            sample_name = cdf_path.stem.removesuffix('_FID1A')
            if allowed_names is not None and sample_name not in allowed_names:
                continue
            metadata = {'AcqMethodName': None, 'InjectionAcqDateTime': None, 'InstrumentName': None,
                        'InjectorPosition': None, 'SampleDescription': None, 'SampleName': sample_name,
                        'SampleType': 'Sample', 'VialNumber': None, 'RawDataFileName': cdf_path.name}
            xy_array = Agilent_FID_Parser.__read_cdf_file(cdf_path)
            injection = FID_Injection(metadata, xy_array, solvent_delay, pos=pos)
            injections[injection.sample_name] = injection

        sequence_metadata = {'sequence_name': Path(raw_directory).stem, 'instrument_name': None, 'instrument': None}
        return sequence_metadata, injections

    @staticmethod
    def __get_xmlns(element:ET.Element) -> str:

        '''Takes in an element, returns the xml-namespace the element or an empty string if no xlmns is found.'''

        m = re.match(r'\{.*\}', element.tag)
        return m.group(0) if m else ''

    @staticmethod
    def __get_root(path: str) -> ET.Element:

        '''Takes in the path to a file in xml format, returns the root node of the xml tree.'''

        tree = ET.parse(path)
        return tree.getroot()

    @staticmethod
    def __get_doc_info(root: ET.Element, xlmns: str) -> ET.Element:

        '''Takes in the root node and xlmns of the xml tree, returns the DocInfo node of the xml tree.'''

        doc = root.find(f'{xlmns}Doc')
        doc_info = doc.find(f'{xlmns}DocInfo')
        return doc_info

    @staticmethod
    def __get_sequence_metadata(root: ET.Element, xlmns: str) -> dict:

        '''
        Returns a dict containing the metadata (sequence_name and instrument_name and instrument) of the sequence.

        Args:
            root (ET.Element): Root node of the xml tree
            xlmns (str): xml-namespace of the root node

        Returns:
            dict: A dictionary containing the sequence_name, instrument_name and instrument
        '''

        doc_info = Agilent_FID_Parser.__get_doc_info(root, xlmns)
        sequence_metadata = {}
        sequence_metadata['sequence_name'] = Agilent_FID_Parser.__get_sequence_name(doc_info, xlmns)
        try:
            name, model = Agilent_FID_Parser.__get_instrument_info(root, xlmns)
        except (AttributeError, TypeError):
            name = Agilent_FID_Parser.__get_instrument_info_fallback(root)
            model = None
        sequence_metadata['instrument_name'] = name
        sequence_metadata['instrument'] = model
        return sequence_metadata

    @staticmethod
    def __get_sequence_name(doc_info:ET.Element, xlmns:str) -> str:

        '''Takes in the DocInfo node of the xml tree and the xlmns and returns the name of the sequence.'''

        sequence_name = doc_info.find(f'{xlmns}Description').text
        return sequence_name

    @staticmethod
    def __get_instrument_info(root:ET.Element, xlmns:str) -> tuple[str,str]:

        '''Takes in the root node of the xml tree and the xlmns and returns the name and model of the instrument.'''

        doc = root.find(f'{xlmns}Doc')
        content = doc.find(f'{xlmns}Content')
        resources = content.find(f'{xlmns}Resources')
        instrument = resources.find(f'{xlmns}Instrument')
        name = instrument.find(f'{xlmns}Name').text
        module = instrument.find(f'{xlmns}Module')
        model = module.find(f'{xlmns}Name').text
        return name, model

    @staticmethod
    def __get_instrument_info_fallback(root: ET.Element) -> Optional[str]:
        """Returns the instrument name read from the first ``InjectionMetaData`` element.

        Used as a fallback when the legacy ``Doc/Content/Resources/Instrument``
        XPath does not exist in the acaml (observed in OpenLab CDS 2.8 SuperGC
        exports).

        Args:
            root: Root node of the acaml xml tree.

        Returns:
            The ``InstrumentName`` attribute of the first ``InjectionMetaData``
            element, or None if no such element is found.
        """
        try:
            metadata_array = list(root.iter('ArrayOfInjectionMetaData'))[0]
            first_injection = list(metadata_array.iter('InjectionMetaData'))[0]
            return first_injection.attrib.get('InstrumentName')
        except (IndexError, KeyError):
            return None


    @staticmethod
    def __get_injections_metadata(root:ET.Element, exclude_blanks:bool = True, sample_filter: Optional[Iterable[str]] = None) -> dict:

        '''
        Returns a dict containing the names of all injections in a sequence and their corresponding metadata.
        The metadata is stored in a dict containing the AcqMethodName, DaMethodName, InjectionAcqDateTime,
        InjectorPosition, SampleDescription, SampleName, SampleType, VialNumber, RawDataFileName.

        Args:
            root (ET.Element): Root node of the xml tree
            exclude_blanks (bool, optional): If True, injections with SampleType 'Blank' are excluded. Defaults to True.
            sample_filter: Optional iterable of allowed ``SampleName`` values. Injections whose ``SampleName`` is not in this set are skipped.
                Defaults to None (no extra filtering).

        Returns:
            dict: A dictionary of all injections of a sequence and their corresponding metadata as a dict
        '''

        allowed_names = set(sample_filter) if sample_filter is not None else None
        injections_metadata = {}
        metadata_array = list(root.iter(f'ArrayOfInjectionMetaData'))[0]
        for metadata_element in metadata_array.iter('InjectionMetaData'):
            metadata = Agilent_FID_Parser.__get_injection_metadata(metadata_element)
            if exclude_blanks:
                if metadata['SampleType'] == 'Blank':
                    continue
            if allowed_names is not None and metadata['SampleName'] not in allowed_names:
                continue

            injections_metadata[metadata['SampleName']] = metadata
        return injections_metadata

    @staticmethod
    def __get_injection_metadata(metadata_element:ET.Element) -> dict:

        '''
        Takes in an array of injection metadata and returns a dict containing the metadata of the injection.

        Args:
            metadata_element (ET.Element): Array of injection metadata

        Returns:
            dict: A dictionary containing the AcqMethodName, DaMethodName, InjectionAcqDateTime,
            InjectorPosition, SampleDescription, SampleName, SampleType, VialNumber of the injection.
        '''

        metadata = {'AcqMethodName': None, 'InjectionAcqDateTime': None, 'InstrumentName': None,
                    'InjectorPosition': None, 'SampleDescription': None, 'SampleName': None, 'SampleType': None,
                    'VialNumber': None, 'RawDataFileName': None}
        for key in metadata.keys():
            if key == 'InjectionAcqDateTime':
                metadata[key] = Agilent_FID_Parser.__convert_to_datetime(metadata_element.attrib[key])
            else:
                metadata[key] = metadata_element.attrib.get(key)
        for child in metadata_element.findall('SampleOrderNumber'):
            metadata['SampleOrderNumber'] = child.attrib['val']
        return metadata

    @staticmethod
    def __convert_to_datetime(inj_acq_time:str) -> datetime:

        '''Takes in the injection aquisition time extracted from the acaml file as a string, returns a corresponding
        datetime object.'''

        date = re.search(r'\d{4}-\d{2}-\d{2}', inj_acq_time).group().split('-')
        time = re.search(r'\d{2}:\d{2}:\d{2}', inj_acq_time).group().split(':')
        date_time = datetime(int(date[0]), int(date[1]), int(date[2]), int(time[0]), int(time[1]), int(time[2]))
        return date_time

    @staticmethod
    def __resolve_file_source(raw_directory: str, file_source: str) -> str:
        """Resolves ``file_source='auto'`` to a concrete strategy.

        Args:
            raw_directory: Path to the result directory.
            file_source: One of ``'csv'``, ``'cdf'`` or ``'auto'``.

        Returns:
            ``'cdf'`` if ``file_source`` is ``'auto'`` and an ``AIA/``
            subdirectory exists, ``'csv'`` if ``file_source`` is ``'auto'`` and
            no ``AIA/`` subdirectory exists, otherwise ``file_source`` unchanged.

        Raises:
            ValueError: If ``file_source`` is not one of the recognized values.
        """
        if file_source not in ('csv', 'cdf', 'auto'):
            raise ValueError(
                f"file_source must be one of 'csv', 'cdf', 'auto'; got "
                f"{file_source!r}"
            )
        if file_source == 'auto':
            return 'cdf' if (Path(raw_directory) / 'AIA').exists() else 'csv'
        return file_source

    @staticmethod
    def __get_xy_arrays(
            xy_directory: str,
            injections_metadata: dict,
            file_source: str = 'csv',
    ) -> dict:
        """Returns a dict containing the xy arrays for all injections in the metadata dict.

        Dispatches to one of two readers based on ``file_source``:
          * ``'csv'``: prefer exact ``<RawDataFileName_stem>.dx_FID1A.CSV``,
            fall back to legacy glob.
          * ``'cdf'``: read ``AIA/<RawDataFileName_stem>_FID1A.cdf`` via
            the ANDI-Chromatography reader.
          * ``'auto'``: resolved to one of the above by ``__resolve_file_source``.

        Args:
            xy_directory: Path to the result directory.
            injections_metadata: Dict of injection metadata keyed by sample
                name (each entry is itself a dict; see
                ``__get_injection_metadata``).
            file_source: Strategy for locating per-injection FID signal files.

        Returns:
            A dict of xy arrays keyed by sample name. Injections whose signal
            file cannot be located are silently skipped (with a warning print).
        """
        resolved_source = Agilent_FID_Parser.__resolve_file_source(xy_directory, file_source)
        xy_arrays = {}
        for sample_name, metadata in injections_metadata.items():
            xy_array = Agilent_FID_Parser.__get_fid_data_for_injection(
                xy_directory, metadata, resolved_source,
            )
            if xy_array is not None:
                xy_arrays[sample_name] = xy_array
        return xy_arrays

    @staticmethod
    def __get_fid_data_for_injection(
            raw_directory: str,
            metadata: dict,
            resolved_source: str,
    ) -> Optional[np.ndarray]:
        """Returns the FID xy array for a single injection.

        Args:
            raw_directory: Path to the result directory.
            metadata: Injection metadata dict (must contain ``SampleName``;
                ``RawDataFileName`` is optional but enables exact-path lookup).
            resolved_source: Either ``'csv'`` or ``'cdf'`` (must already be
                resolved; ``'auto'`` is not accepted here).

        Returns:
            A 2xN numpy array ``[time_min, intensity]``, or None if the file
            cannot be located.
        """
        raw_path = Path(raw_directory)
        sample_name = metadata['SampleName']
        raw_data_name = metadata.get('RawDataFileName') or sample_name
        raw_data_stem = Path(raw_data_name).stem

        if resolved_source == 'cdf':
            cdf_path = raw_path / 'AIA' / f'{raw_data_stem}_FID1A.cdf'
            if not cdf_path.exists():
                print(f'Warning: No matching FID CDF for sample {sample_name} at {cdf_path}')
                return None
            return Agilent_FID_Parser.__read_cdf_file(cdf_path)

        # Legacy CSV/.xy path. Prefer exact RawDataFileName-based lookup; fall
        # back to the original glob-based behavior so callers without
        # ``RawDataFileName`` in their acaml still work.
        exact_csv = raw_path / f'{raw_data_stem}.dx_FID1A.CSV'
        exact_xy = raw_path / f'{raw_data_stem}.xy'
        if exact_csv.exists():
            return FID_Base_Parser.read_xy_array(exact_csv)
        if exact_xy.exists():
            return FID_Base_Parser.read_xy_array(exact_xy)
        # Legacy glob fallback (preserves prior behavior for older datasets).
        supported_formats = {'.xy', '.CSV', '.csv', '.XY'}
        file_paths = [
            fp for fp in raw_path.glob(f'*{sample_name}*')
            if fp.suffix in supported_formats
        ]
        if not file_paths:
            print(f'Warning: No matching file for sample {sample_name}')
            return None
        return FID_Base_Parser.read_xy_array(file_paths[0])

    @staticmethod
    def __read_cdf_file(cdf_path: Path) -> np.ndarray:
        """Reads (.cdf) file and returns a 2xN xy array.

        The time axis is reconstructed from the file's
        ``actual_sampling_interval`` (and optional ``actual_delay_time``) and
        converted to minutes to match the convention of the legacy .xy/.CSV
        readers.

        Args:
            cdf_path: Path to the .cdf file.

        Returns:
            A 2xN numpy array ``[time_min, intensity]``, matching the shape
            returned by ``FID_Base_Parser.read_xy_array``.
        """
        dataset = nc.Dataset(cdf_path, 'r')
        try:
            intensities = np.asarray(dataset.variables['ordinate_values'][:],
                                     dtype=float)
            sampling_interval_s = float(np.asarray(
                dataset.variables['actual_sampling_interval'][:]).item())
            try:
                delay_s = float(np.asarray(
                    dataset.variables['actual_delay_time'][:]).item())
            except KeyError:
                delay_s = 0.0
        finally:
            dataset.close()

        n_points = len(intensities)
        times_s = delay_s + np.arange(n_points) * sampling_interval_s
        times_min = times_s / 60.0
        return np.vstack([times_min, intensities])

    @staticmethod
    def __initialize_injections(injections_metadata: dict, xy_arrays: dict, solvent_delay:float, pos:bool=False) -> dict[str:FID_Injection]:

        '''Takes in a dict of injection metadata and xy_arrays, returns a dict of injection objects accessed by the
        sample name.'''

        injections = {}
        for name, metadata in injections_metadata.items():
            injection = FID_Injection(metadata, xy_arrays[name], solvent_delay, pos=pos)
            injections[injection.sample_name] = injection
        return injections



if __name__ == '__main__':
    pass
