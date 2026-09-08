import numpy as np

class Analysis_Settings:

    '''
    Analysis_Settings class for storing information about parameters for raw data processing.

    Attributes:
        sn (int): Signal to noise ratio for peak detection.
        time_range (tuple): Time range for peak detection.
        width (int, float): Width for peak detection.
        prominence_ms (int): Prominence for peak detection via MS.
        prominence_fid (int): Prominence for peak detection via FID.
        trace_prominence (int): Prominence for peak detection on individual mass traces.
        height (int): Height for peak detection.
        savgol_window (int): Window size for Savitzky-Golay filter.
        max_half_window (int): Maximum half window size for Snip algorithm.
        boarder_threshold (int): Threshold for boarder detection.
        boarder_window (int): Window size for boarder detection.
        max_isotopic_diff (float): Max deviation of isotopic peak intensity for analyte assignment.
        min_rel_intensity (float): Minimum relative intensity for m/z trace to be considered for analyte assignment.
        min_mz_fraction (float): Minimum fraction of a spectrum's maximum m/z that the parent ion must exceed to
        be considered for analyte assignment.
        scan_rate (float): Scan rate of chromatogram.
        _resolved (dict): Values handed out by pop since the recording decorator last cleared it,
        keyed by setting name. Not a setting itself: it is rejected by update and by pop.
    '''

    sn: int
    time_range: tuple|None
    width: int|float|None
    prominence_ms: int | None
    prominence_fid: int | None
    trace_prominence: int|None
    height: int|None
    savgol_window: int|None
    max_half_window: int|None
    boarder_threshold: int|None
    boarder_window: int|None
    max_isotopic_diff: float|None
    min_rel_intensity: float|None
    min_mz_fraction: float|None
    scan_rate: float
    _resolved: dict


    __slots__ = 'sn', 'time_range', 'width', 'prominence_ms', 'prominence_fid', 'trace_prominence', 'height', \
                'savgol_window', 'max_half_window', 'boarder_threshold', 'boarder_window', 'max_isotopic_diff', 'min_rel_intensity', 'min_mz_fraction', 'scan_rate', '_resolved'

    def __init__(self, chromatogram:np.ndarray):
        self.sn = 5
        self.time_range = None
        self.width = None
        self.prominence_ms = None
        self.prominence_fid = None
        self.trace_prominence = None
        self.height = None
        self.savgol_window = None
        self.max_half_window = None
        self.boarder_threshold = None
        self.boarder_window = None
        self.max_isotopic_diff = None
        self.min_rel_intensity = None
        self.min_mz_fraction = None
        self.scan_rate = chromatogram[0, 2] - chromatogram[0, 1]
        self._resolved = {}

    def __str__(self) -> str:
        return f'Analysis_Settings:\nSignal to Noise Ratio: {self.sn}\nTime Range: {self.time_range}\n' \
               f'Width: {self.width}\nProminence MS: {self.prominence_ms}\nProminence MS: {self.prominence_ms}\n' \
               f'Trace Prominence: {self.trace_prominence}\nHeight: {self.height}\n' \
               f'Savitzky-Golay Window: {self.savgol_window}\nMax Half Window: {self.max_half_window}\n' \
               f'Boarder Threshold: {self.boarder_threshold}\nBoarder Window: {self.boarder_window}\n' \
               f'Max Isotopic Diff: {self.max_isotopic_diff}\nMin Relative Intensity: {self.min_rel_intensity}\n' \
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

    def pop(self, key:str, default):

        '''
        Returns the configured value for a setting if one is set and the given default otherwise,
        recording the returned value in _resolved.

        Note this does not remove anything; the name is historic. Most thresholds in the library are
        never configured and are computed from the signal at the call site, so the value returned
        here is the only record of what an algorithm actually used.

        Args:
            key (str): Name of the setting.
            default: Value to use when the setting is not configured.

        Returns:
            The value the caller will use for the setting.
        '''

        if key not in self.__slots__ or key.startswith('_'):
            raise KeyError(f'"{key}" is not a valid setting.')
        value = getattr(self, key)
        if not value:
            value = default
        self._resolved[key] = value
        return value

    def __setstate__(self, state:tuple) -> None:

        '''
        Restores Analysis_Settings from its pickled state, defaulting attributes added since the
        file was written and dropping ones removed since.

        Settings pickled before pop recorded resolved values carry no _resolved entry, and every pop
        on them would raise AttributeError without the default; ones pickled while indices_range was
        still a slot carry a value for it, and setting it would raise AttributeError without the
        skip. The settings object is nested inside a pickled injection and restores itself, so
        Injection.__setstate__ cannot cover either case.
        '''

        _, slots = state
        for name, value in (slots or {}).items():
            # Skip attributes that are no longer slots: a .pkl is coupled to the class layout
            # (architecture 8), and indices_range was removed once the window came to be derived
            # from the chromatogram at the call site. hasattr on the type finds slot descriptors
            # across the MRO, which self.__slots__ does not.
            if hasattr(type(self), name):
                setattr(self, name, value)
        if not hasattr(self, '_resolved'):
            self._resolved = {}


    def __check_settings(self, setting:str, value:int|float|tuple) -> bool:

        '''
        Returns True if setting is valid.

        Args:
            setting (str): Setting to be checked.
            value (int, float, tuple): Value to be checked.

        Returns:
            bool: True if setting is valid.
        '''

        options = {'sn': int, 'time_range': tuple, 'width': (int, float), 'prominence_ms': (int, float), 'prominence_fid': int,
                   'trace_prominence': int, 'height': int, 'savgol_window': int, 'max_half_window': int,
                   'boarder_threshold': int, 'boarder_window': int, 'max_isotopic_diff': float, 'min_rel_intensity': float,
                   'min_mz_fraction': float}
        if setting in options.keys():
            if isinstance(value, options[setting]):
                return True
            else:
                raise TypeError(f'"{setting}" is expected to be {type(options[setting])} not {type(value)}.')
        else:
            raise KeyError(f'"{setting}" is not a valid setting.')
