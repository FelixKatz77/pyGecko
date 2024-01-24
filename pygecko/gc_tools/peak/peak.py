import numpy as np

from pygecko.gc_tools.analyte import Analyte


class Peak:
    '''
    Base class for all Peaks.

    Attributes:
        rt (float): Retention time of the peak.
        ri (None|float): Retention index of the peak.
        area (None|float): Area of the peak. Defaults to None.
        height (None|float): Height of the peak.
        width (None|float): Width of the peak.
        boarders (None|np.ndarray): Boarders of the peak.
        flag (None|str): Flag of the peak.
        analyte (None|Analyte): Analyte of the peak.
    '''

    rt: float
    ri: None|float
    area: None|float
    height: None|float
    width: None|float
    boarders: None|np.ndarray
    flag: None|str
    analyte: None|Analyte

    __slots__ = 'rt', 'ri', 'area', 'height', 'width', 'boarders', 'flag', 'analyte'

    def __init__(self, rt: float, height: float, width: float, borders: np.ndarray, area=None):
        self.rt = rt
        self.ri = None
        self.area = area
        self.height = height
        self.width = width
        self.boarders = borders
        self.flag = None
        self.analyte = None



