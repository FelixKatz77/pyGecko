import logging
from pathlib import Path

import numpy as np

from pygecko.gc_tools import Chromatogram, FID_Injection, FID_Sequence, RI_Calibration

logger = logging.getLogger(__name__)


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

        logger.info('Loading GC-FID sequence...')
        xy_directory = Path(xy_directory)
        # Matched on the lower-cased suffix of each entry, so every case variant is found on a
        # case-sensitive filesystem and each file is listed once on a case-insensitive one.
        supported_formats = {'.xy', '.csv'}
        xy_files = sorted(entry for entry in xy_directory.iterdir()
                          if entry.suffix.lower() in supported_formats)
        injections = {}
        for xy_file in xy_files:
            injection = FID_Base_Parser.load_injection(xy_file, solvent_delay, pos=pos)
            injections[injection.sample_name] = injection
        logger.info('Sequence loaded with %d injections.', len(injections))
        return FID_Sequence({}, injections)

    @staticmethod
    def load_ri_calibration(xy_file: Path|str, solvent_delay, c_count: int, rt: float) -> RI_Calibration:

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
        injection = FID_Injection({'SampleName': sample_name}, Chromatogram(*xy_array, kind='FID'), solvent_delay)
        injection.record_step('FID_Base_Parser.load_injection',
                              {'xy_file': str(xy_file), 'solvent_delay': solvent_delay})
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
        injection = FID_Injection({'SampleName':sample_name}, Chromatogram(*xy_array, kind='FID'), solvent_delay,
                                  pos=pos)
        injection.record_step('FID_Base_Parser.load_injection',
                              {'xy_file': str(xy_file), 'solvent_delay': solvent_delay, 'pos': pos})
        return injection


    @staticmethod
    def read_xy_array(path:Path) -> np.ndarray:

        '''Takes in the path to a tab-separated .xy or comma-separated .csv file, in any case, returns the
        xy_array.'''

        delimiter = {'.xy': '\t', '.csv': ','}.get(path.suffix.lower())
        if delimiter is None:
            raise ValueError(f'Cannot read {path.name}: unsupported suffix {path.suffix!r}, expected .xy or .csv.')
        return np.transpose(np.loadtxt(path, delimiter=delimiter, converters={0: float}))


