import numpy as np

from pygecko.gc_tools.peak import Peak


class FID_Peak(Peak):
    '''
    FID_Peak class for storing information about peaks in FID chromatograms.

    Attributes:
        rt (float): Retention time of the peak.
        ri (None|float): Retention index of the peak.
        area (None|float): Area of the peak. Defaults to None.
        flags (list[Any|str]): List of flags of the peak.
        analyte (None|Analyte): Analyte of the peak.
        height (float): Height of the peak.
        boarders (np.ndarray): Boarders of the peak.
        width (float): Width of the peak.
        boarders (np.ndarray): Boarders of the peak.
    '''

    def __init__(self, rt: float, height: float, width: float, boarders: np.ndarray, area: float, flag: str|None):
        super().__init__(rt, height, width, boarders, area=area)
        if flag:
            self.flags = [flag]
