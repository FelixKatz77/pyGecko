import numpy as np
from pathlib import Path

from pygecko.gc_tools import FID_Injection, FID_Sequence, RI_Calibration


class FID_Base_Parser:

    @staticmethod
    def load_sequence(xy_directory: Path|str, solvent_delay:float, pos:bool=False) -> FID_Sequence:

        '''
        Returns an FID_Sequence object.

        Args:
            xy_directory (Path|str): Path to a directory containing xy-files.
            solvent_delay (float): Solvent delay of the injections.

        Returns:
            FID_Sequence: An FID_Sequence object.
        '''

        xy_directory = Path(xy_directory)
        supported_formats = ['.xy', '.CSV']
        xy_files = []
        for file_format in supported_formats:
            xy_files.extend(xy_directory.glob(f'*{file_format}'))
        injections = {}
        for xy_file in xy_files:
            injection = FID_Base_Parser.load_injection(xy_file, solvent_delay, pos=pos)
            injections[injection.sample_name] = injection
        return FID_Sequence({}, injections)

    @classmethod
    def load_ri_calibration(cls, xy_file: Path|str, solvent_delay, c_count: int, rt: float) -> RI_Calibration:

        '''
        Returns an RI_Calibration object.

        Args:
            xy_file (Path|str): Path to a xy_file.
            solvent_delay (float): Solvent delay of the injection.
            c_count (int): Carbon count of the alkane the retention time is provided for.
            rt (float): Retention time of the alkane the c_count is provided for.

        Returns:
            RI_Calibration: An RI_Calibration object.
        '''

        xy_file = Path(xy_file)
        xy_array = FID_Base_Parser.read_xy_array(xy_file)
        sample_name = xy_file.stem.split('.')[0]
        injection = FID_Injection({'SampleName': sample_name}, xy_array, solvent_delay)
        return RI_Calibration(injection, c_count, rt)

    @staticmethod
    def load_injection(xy_file: Path|str, solvent_delay:float, pos:bool=False) -> FID_Injection:
        '''
        Returns an FID_Injection object.

        Args:
            xy_file (Path|str): Path to a xy_file.
            solvent_delay (float): Solvent delay of the injection.

        Returns:
            FID_Injection: An FID_Injection object.
        '''

        xy_file = Path(xy_file)
        xy_array = FID_Base_Parser.read_xy_array(xy_file)
        sample_name = xy_file.stem.split('.')[0]
        injection = FID_Injection({'SampleName':sample_name}, xy_array, solvent_delay, pos=pos)
        return injection


    @staticmethod
    def read_xy_array(path:Path) -> np.ndarray|None:

        '''Takes in the path to a xy-file, returns the xy_array.'''

        if path.suffix == '.xy':
            array = np.transpose(np.loadtxt(path, delimiter='\t'))
            return array
        elif path.suffix == '.CSV':
            array = np.transpose(np.loadtxt(path, delimiter=',', converters={0: float}))
            return array
        else:
            print(f'Cannot read {path.name}: File format not supported.')
            return None


