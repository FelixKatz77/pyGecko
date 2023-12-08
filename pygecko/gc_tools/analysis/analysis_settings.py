import numpy as np
from gc_tools.utilities import Utilities

class Analysis_Settings:

    '''
    Analysis_Settings class for storing information about parameters for raw data processing.

    Attributes:
        sn (int): Signal to noise ratio for peak detection.
        time_range (tuple): Time range for peak detection.
        indices_range (list): Indices range for peak detection.
        width (int, float): Width for peak detection.
        prominence (int): Prominence for peak detection.
        trace_prominence (int): Prominence for peak detection on individual mass traces.
        height (int): Height for peak detection.
        savgol_window (int): Window size for Savitzky-Golay filter.
        max_half_window (int): Maximum half window size for Snip algorithm.
        boarder_threshold (int): Threshold for boarder detection.
        boarder_window (int): Window size for boarder detection.
        scan_rate (float): Scan rate of chromatogram.
    '''

    sn: int
    time_range: tuple|None
    indices_range: list[int|None]
    width: int|float|None
    prominence: int|None
    trace_prominence: int|None
    height: int|None
    savgol_window: int|None
    max_half_window: int|None
    boarder_threshold: int|None
    boarder_window: int|None
    scan_rate: float

    __slots__ = 'sn', 'time_range', 'indices_range', 'width', 'prominence', 'trace_prominence', 'height', \
                'savgol_window', 'max_half_window', 'boarder_threshold', 'boarder_window', 'scan_rate'

    def __init__(self, chromatogram:np.ndarray):
        self.sn = 5
        self.time_range = None
        self.indices_range = self.__set_indices_range()
        self.width = None
        self.prominence = None
        self.trace_prominence = None
        self.height = None
        self.savgol_window = None
        self.max_half_window = None
        self.boarder_threshold = None
        self.boarder_window = None
        self.scan_rate = chromatogram[0, 2] - chromatogram[0, 1]

    def __str__(self) -> str:
        return f'Analysis_Settings:\nSignal to Noise Ratio: {self.sn}\nTime Range: {self.time_range}\nWidth: {self.width}\n' \
               f'Prominence: {self.prominence}\nTrace Prominence: {self.trace_prominence}\nHeight: {self.height}\n' \
               f'Savitzky-Golay Window: {self.savgol_window}\nMax Half Window: {self.max_half_window}\n' \
               f'Boarder Threshold: {self.boarder_threshold}\nBoarder Window: {self.boarder_window}\n' \
               f'Scan Rate: {self.scan_rate}'


    def update(self, **kwargs):

        '''
        Updates settings by setting the keyword arguments after checking.

        Args:
            **kwargs: Keyword arguments for settings.
        '''

        for key, value in kwargs.items():
            if self.__check_settings(key, value):
                setattr(self, key, value)
        self.indices_range = self.__set_indices_range()

    def pop(self, key:str, default):
        if key in self.__slots__:
            value = getattr(self, key)
            if value:
                return value
            else:
                return default
        else:
            raise KeyError(f'"{key}" is not a valid setting.')


    def __check_settings(self, setting:str, value:int|float|tuple) -> bool:

        '''
        Returns True if setting is valid.

        Args:
            setting (str): Setting to be checked.
            value (int, float, tuple): Value to be checked.

        Returns:
            bool: True if setting is valid.
        '''

        options = {'sn': int, 'time_range': tuple, 'width': (int, float), 'prominence': int, 'trace_prominence': int,
                   'height': int, 'savgol_window': int, 'max_half_window': int, 'boarder_threshold': int,
                   'boarder_window': int}
        if setting in options.keys():
            if isinstance(value, options[setting]):
                return True
            else:
                raise TypeError(f'"{setting}" is expected to be {type(options[setting])} not {type(value)}.')
        else:
            raise KeyError(f'"{setting}" is not a valid setting.')

    def __set_indices_range(self) -> list[int|None]:

        '''
        Returns indices range for peak detection.

        Returns:
            list: Indices range for peak detection.
        '''

        if not self.time_range:
            indices_range = [None, None]
        else:
            indices_range = Utilities.convert_time_to_scan(self.time_range, self.scan_rate)
        if indices_range[0] is None:
            indices_range[0] = 0
        if indices_range[1] is None:
            indices_range[1] = None
        return indices_range