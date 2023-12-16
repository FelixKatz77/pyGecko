from pygecko.gc_tools.analyte import Analyte


class Peak:
    '''
    Base class for all Peaks.

    Attributes:
        rt (float): Retention time of the peak.
        ri (None|float): Retention index of the peak.
        area (None|float): Area of the peak. Defaults to None.
        height (float): Height of the peak.
        flag (None|str): Flag of the peak.
        analyte (None|Analyte): Analyte of the peak.
    '''

    rt: float
    ri: None|float
    area: None|float
    flag: None|str
    analyte: None|Analyte

    __slots__ = 'rt', 'ri', 'area', 'height', 'flag', 'analyte'

    def __init__(self, rt: float, height: float, area=None):
        self.rt = rt
        self.ri = None
        self.area = area
        self.height = height
        self.flag = None
        self.analyte = None



