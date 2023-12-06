import numpy as np
from pathlib import Path

from gc_tools import FID_Injection, FID_Sequence


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
        xy_files = xy_directory.glob('*.xy')
        injections = {}
        for xy_file in xy_files:
            injection = FID_Base_Parser.load_injection(xy_file, solvent_delay, pos=pos)
            injections[injection.sample_name] = injection
        return FID_Sequence({}, injections)


    @staticmethod
    def load_injection(xy_file: Path, solvent_delay:float, pos:bool=False) -> FID_Injection:
        '''
        Returns an FID_Injection object.

        Args:
            xy_file (Path): Path to a xy_file.
            solvent_delay (float): Solvent delay of the injection.

        Returns:
            FID_Injection: An FID_Injection object.
        '''

        xy_array = FID_Base_Parser.__read_xy_array(Path(xy_file))
        sample_name = xy_file.stem
        injection = FID_Injection({'SampleName':sample_name}, xy_array, solvent_delay, pos=pos)
        return injection


    @staticmethod
    def __read_xy_array(path:Path) -> np.ndarray:

        '''Takes in the path to a xy-file, returns the xy_array.'''

        array = np.loadtxt(path)
        array = np.transpose(array)
        return array