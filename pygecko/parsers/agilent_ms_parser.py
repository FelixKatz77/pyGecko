import tempfile
from pathlib import Path
from xml.etree import ElementTree as ET

from pygecko.gc_tools import MS_Injection, MS_Sequence, RI_Calibration


from pygecko.parsers.ms_base_parser import MS_Base_Parser
from .utilities import find_directories_with_extension


class Agilent_MS_Parser:

    '''
    A class wrapping functions to parse raw data from an Agilent GC/MS into MS_Sequence objects.
    '''

    @classmethod
    def load_sequence(cls, sequence_directory: str, pos:bool=False) -> MS_Sequence:

        '''
        Returns an MS_Sequence object.

        Args:
            sequence_directory (str): Path to the directory containing the raw data.
            chrom_path (str): Path to the directory containing the chromatograms.
            pos (bool, optional): Indicates if plate position is given in the injection names. Defaults to False.

        Returns:
            MS_Sequence: An MS_Sequence object
        '''

        sequence_metadata, injections = cls.__load_sequence_data(sequence_directory, pos=pos)
        return MS_Sequence(sequence_metadata, injections)

    @classmethod
    def load_ri_calibration(cls, raw_directory: str, c_count:int, rt:float) -> RI_Calibration:

        '''
        Returns an RI_Calibration object.

        Args:
            raw_directory (str): Path to the directory containing the raw data.
            c_count (int): Number of carbon atoms for as specific alkane present in the standard.
            rt (float): Retention time of the alkane the c_count is provided for.

        Returns:
            RI_Calibration: An RI_Calibration object
        '''

        injection = cls.__load_injection_data(raw_directory)
        return RI_Calibration(injection, c_count, rt)

    @classmethod
    def load_injection(cls, raw_directory: str) -> MS_Injection:

        '''
        Returns an MS_Injection object.

        Args:
            raw_directory (str): Path to the directory containing the raw data.

        Returns:
            MS_Injection: An MS_Injection object
        '''

        injection = cls.__load_injection_data(raw_directory)
        return injection

    @staticmethod
    def __load_injection_data(raw_directory:str, temp_dir:tempfile.TemporaryDirectory=None, pos:bool=False) -> MS_Injection:

        '''
        Returns an MS_Injection object.

        Args:
            raw_directory (str): Path to the directory containing the raw data.
            chrom_path (str): Path to the directory containing the chromatogram.

        Returns:
            MS_Injection: An MS_Injection object
        '''
        xml_file = Agilent_MS_Parser.__get_sequence_xml_file(raw_directory)
        injection_metadata, sequence_info = Agilent_MS_Parser.__get_injection_metadata(xml_file)
        injection = Agilent_MS_Parser.__initialize_injection(Path(raw_directory), injection_metadata, pos=pos, temp_dir=temp_dir)
        return injection
    @staticmethod
    def __load_sequence_data(sequence_directory: str, pos:bool=False) -> (dict[str:str], dict[str:MS_Injection]):

        '''
        Returns a dict containing the sequence metadata and a dict containing the injections.

        Args:
            sequence_directory (str): Path to the directory containing the raw data.

        Returns:
            tuple[dict, dict[str:MS_Injection]]: A tuple containing a dict containing the sequence metadata and a
            dict containing the injections.
        '''

        temp_dir = tempfile.TemporaryDirectory()
        raw_directories = find_directories_with_extension(sequence_directory, '.D')
        injections = {}
        for raw_directory in raw_directories:
            injection = Agilent_MS_Parser.__load_injection_data(raw_directory, temp_dir=temp_dir, pos=pos)
            injections[injection.sample_name] = injection
        sequence_metadata = Agilent_MS_Parser.__get_sequence_metadata(sequence_directory, injection.instrument_name)
        return sequence_metadata, injections

    @staticmethod
    def __get_sequence_metadata(raw_directory:str, instrument_name:str) -> dict[str:str]:

        '''
        Takes in the path to the directory containing the raw data and the name of the instrument, returns a dict
        containing the sequence metadata.
        '''

        sequence_metadata = {}
        sequence_metadata['sequence_name'] = Path(raw_directory).name
        sequence_metadata['instrument_name'] = instrument_name
        return sequence_metadata

    @staticmethod
    def __get_injections_metadata(raw_directory: str) -> (dict[str:dict], str):

        '''
        Takes in the path to the directory containing the raw data, returns a tuple containing a dict with each
        injection's metadata (AcqMethodFileName, SampleInformation, SampleName, SampleType, Vial) accessed by the
        sample name and the name of the instrument.
        '''

        sequence_xml_files = Agilent_MS_Parser.__get_sequence_xml_file(raw_directory)
        injections_metadata = {}
        for xml_file in sequence_xml_files:
            metadata, sequence_info = Agilent_MS_Parser.__get_injection_metadata(xml_file)
            injections_metadata[metadata['SampleName']] = metadata
        instrumnet_name = Agilent_MS_Parser.__get_instrument_info(sequence_info)
        return injections_metadata, instrumnet_name

    @staticmethod
    def __get_injection_metadata(xml_file:str) -> dict[str:str]:

        '''
        Takes in the sequence.xml file, returns a dict
        containing the injection metadata.
        '''


        root = ET.parse(xml_file).getroot()
        sequence_info = root.find(f'Sequence')
        metadata = {'AcqMethodFileName': None, 'InstrumentName': None, 'SampleInformation': None, 'SampleName': None, 'SampleType': None,
                    'Vial': None}
        for key in metadata:
            metadata[key] = sequence_info.find(key).text
        metadata['AcqMethodName'] = metadata.pop('AcqMethodFileName')
        metadata['SampleDescription'] = metadata.pop('SampleInformation')
        metadata['VialNumber'] = metadata.pop('Vial')
        return metadata, sequence_info


    @staticmethod
    def __get_instrument_info(sequence_info: ET.Element) -> str:

        '''
        Takes in the Sequence node of the xml tree, returns the name of the instrument.
        '''

        instrument_name = sequence_info.find('InstrumentName').text
        return instrument_name

    @staticmethod
    def __get_sequence_xml_file(directory: str) -> str:

        '''
        Takes in the path to the directory containing the raw data, returns a list containing the paths to all
        sequence.xml files as strings.
        '''

        xml_file =  str(list(Path(directory).rglob('sequence.xml'))[0])
        return xml_file

    @staticmethod
    def __initialize_injection(path:Path, metadata:dict[str, str], pos:bool=False, temp_dir:tempfile.TemporaryDirectory=None) -> MS_Injection:

        '''
        Returns an MS_Injection object.

        Args:
            path (Path): Path to the Agilent .D directory.
            metadata (dict[str, str]): Injection metadata.

        Returns:
            MS_Injection: An MS_Injection object
        '''

        scans = MS_Base_Parser.extract_scans_from_raw_data(path, temp_dir=temp_dir)
        chromatogram = MS_Base_Parser.construct_chromatogram_from_scans(scans)
        injection = MS_Injection(metadata, chromatogram, None, scans, pos=pos)
        return injection






