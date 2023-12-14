import xml.etree.ElementTree as ET
import re
import numpy as np
from pathlib import Path
from datetime import datetime

from pygecko.parsers.fid_base_parser import FID_Base_Parser
from pygecko.gc_tools import FID_Sequence, FID_Injection, RI_Calibration


class Agilent_FID_Parser:

    '''
    A class wrapping functions to parse raw data from an Agilent GC/FID into FID_Sequence objects.
    '''

    @classmethod
    def load_sequence(cls, raw_directory:str, solvent_delay:float|int, pos:bool=False) -> FID_Sequence:
        '''
        Returns an FID_Sequence object.

        Args:
            raw_directory (str): Path to the directory containing the corresponding xy-files and the acaml-file.
            solvent_delay (float): Retention time of the solvent peak in minutes.
            pos (bool, optional): Indicates if plate position is given in the injection names. Defaults to False.

        Returns:
            FID_Sequence: An FID_Sequence object
        '''

        acaml_file = list(Path(raw_directory).glob('*.acaml'))[0]
        sequence_metadata, injections = cls.__load_sequence_data(acaml_file, raw_directory, solvent_delay, pos=pos)
        return FID_Sequence(sequence_metadata, injections)

    @classmethod
    def load_ri_calibration(cls, raw_directory:str, solvent_delay, c_count:int, rt:float) -> RI_Calibration:

        '''
        Returns an RI_Calibration object.

        Args:
            raw_directory (str): Path to the directory containing the corresponding xy-files and the acaml-file.
            solvent_delay (float): Retention time of the solvent peak in minutes.
            c_count (int): Carbon count of the alkane the retention time is provided for.
            rt (float): Retention time of the alkane the c_count is provided for.

        Returns:
            RI_Calibration: An RI_Calibration object.
        '''

        acaml_file = list(Path(raw_directory).glob('*.acaml'))[0]
        sequence_metadata, injections = cls.__load_sequence_data(acaml_file, raw_directory, solvent_delay)
        return RI_Calibration(sequence_metadata, injections, c_count, rt)

    @classmethod
    def load_injection(cls, raw_path:str, solvent_delay:float|int) -> FID_Injection:

        '''
        Returns an FID_Injection object.

        Args:
            raw_path (str): Path to the directory containing the corresponding xy-file and the acaml-file.
            solvent_delay (float): Retention time of the solvent peak in minutes.

        Returns:
            FID_Injection: An FID_Injection object.
        '''

        acaml_file = list(Path(raw_path).glob('*.acaml'))[0]
        injection = cls.__load_injection_data(acaml_file, raw_path, solvent_delay)
        return injection

    @staticmethod
    def __load_injection_data(acaml_file:Path, xy_file:str, solvent_delay:float|int) -> FID_Injection:

        '''
        Returns an FID_Injection object.

        Args:
            acaml_file (Path): Path to an acaml-file from an Agilent GC.
            xy_file (str): Path to the corresponding xy-file.
            solvent_delay (float): Retention time  of the solvent peak in minutes.

        Returns:
            FID_Injection: An FID_Injection object.
        '''

        root = Agilent_FID_Parser.__get_root(acaml_file)
        metadata_array = list(root.iter(f'ArrayOfInjectionMetaData'))[0]
        metadata_element = list(metadata_array.iter('InjectionMetaData'))[0]
        injection_metadata = Agilent_FID_Parser.__get_injection_metadata(metadata_element)
        xy_array = FID_Base_Parser.read_xy_array(Path(xy_file))
        injection = FID_Injection(injection_metadata, xy_array, solvent_delay)
        return injection


    @staticmethod
    def __load_sequence_data(acaml_file:Path, xy_directory:str, solvent_delay:float|int, pos:bool=False) -> (dict, dict[str:FID_Injection]):

        '''
        Returns a dict containing the sequence metadata and a dict containing the injections.

        Args:
            acaml_file (Path): Path to an acaml-file from an Agilent GC
            xy_directory (str): Path to the directory containing the corresponding xy-files.
            solvent_delay (float): Retention time  of the solvent peak in minutes.
            pos (bool, optional): Indicates if plate position is given in the injection names. Defaults to False.

        Returns:
            tuple[dict, dict[str:FID_Injection]]: A tuple containing a dict containing the sequence metadata and
            a dict containing the injections.
        '''

        root = Agilent_FID_Parser.__get_root(acaml_file)
        xlmns = Agilent_FID_Parser.__get_xmlns(root)
        sequence_metadata = Agilent_FID_Parser.__get_sequence_metadata(root, xlmns)
        injections_metadata = Agilent_FID_Parser.__get_injections_metadata(root)
        xy_arrays = Agilent_FID_Parser.__get_xy_arrays(xy_directory, injections_metadata.keys())
        injections = Agilent_FID_Parser.__initialize_injections(injections_metadata, xy_arrays, solvent_delay, pos=pos)
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
        sequence_metadata['instrument_name'], sequence_metadata['instrument'] = Agilent_FID_Parser.__get_instrument_info(
            root, xlmns)
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
    def __get_injections_metadata(root:ET.Element, exclude_blanks:bool = True) -> dict:

        '''
        Returns a dict containing the names of all injections in a sequence and their corresponding metadata.
        The metadata is stored in a dict containing the AcqMethodName, DaMethodName, InjectionAcqDateTime,
        InjectorPosition, SampleDescription, SampleName, SampleType, VialNumber.

        Args:
            root (ET.Element): Root node of the xml tree
            exclude_blanks (bool, optional): If True, injections with SampleType 'Blank' are excluded. Defaults to True.

        Returns:
            dict: A dictionary of all injections of a sequence and their corresponding metadata as a dict
        '''

        injections_metadata = {}
        metadata_array = list(root.iter(f'ArrayOfInjectionMetaData'))[0]
        for metadata_element in metadata_array.iter('InjectionMetaData'):
            metadata = Agilent_FID_Parser.__get_injection_metadata(metadata_element)
            if exclude_blanks:
                if metadata['SampleType'] == 'Blank':
                    continue
                else:
                    pass
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
                    'VialNumber': None}
        for key in metadata.keys():
            if key == 'InjectionAcqDateTime':
                metadata[key] = Agilent_FID_Parser.__convert_to_datetime(metadata_element.attrib[key])
            else:
                metadata[key] = metadata_element.attrib[key]
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
    def __get_xy_arrays(xy_directory:str, sample_names:iter) -> dict:

        '''Takes in the directory containing the xy-files for injections and the injections' names, returns a dict
        containing the xy_arrays for the injections.'''

        xy_arrays = {}
        path = Path(xy_directory)
        supported_formats = ['.xy', '.CSV']
        for sample_name in sample_names:
            file_paths = list(path.glob(f'{sample_name}*.*'))
            file_paths = list(filter(lambda x: x.suffix in supported_formats, file_paths))
            xy_array = FID_Base_Parser.read_xy_array(file_paths[0])
            if xy_array is not None:
                xy_arrays[sample_name] = xy_array
            else:
                pass
        return xy_arrays


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