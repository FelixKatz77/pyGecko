import numpy as np

from pygecko.gc_tools.peak import Peak
from visualization import Visualization


class MS_Peak(Peak):
    '''
    MS_Peak class for storing information about peaks in MS chromatograms.

    Attributes:
        rt (float): Retention time of the peak.
        ri (Union[None, float]): Retention index of the peak.
        area (Union[None, float]): Area of the peak. Defaults to None.
        flag (Union[None, str]): Flag of the peak.
        analyte (Union[None, Analyte]): Analyte of the peak.
        mass_spectrum (np.ndarray): Mass spectrum of the peak.
    '''

    mass_spectrum: np.ndarray

    __slots__ = 'mass_spectrum'

    def __init__(self, rt: float, mass_spectrum: np.ndarray, area=None):
        super().__init__(rt, area=area)
        self.mass_spectrum = mass_spectrum

    def view_mass_spectrum(self, path:str|None=None, **kwargs) -> None:
        '''
        Plots the mass spectrum of the peak.
        Args:
            **kwargs: Keyword arguments for the visualization.
        '''

        Visualization.view_mass_spectrum(self, path=path, **kwargs)

    def __contains__(self, item):

        '''
        Returns True if the given m/z value is in the peak's mass spectrum and False otherwise.
        '''

        return item in self.mass_spectrum['mz']

