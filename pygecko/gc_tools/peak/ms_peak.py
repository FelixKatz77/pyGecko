import numpy as np
import pandas as pd

from pygecko.gc_tools.peak import Peak
from pygecko.visualization import Visualization


class MS_Peak(Peak):
    '''
    MS_Peak class for storing information about peaks in MS chromatograms.

    Attributes:
        rt (float): Retention time of the peak.
        ri (Union[None, float]): Retention index of the peak.
        area (Union[None, float]): Area of the peak. Defaults to None.
        height (float): Height of the peak.
        width (float): Width of the peak.
        boarders (np.ndarray): Boarders of the peak.
        flags (list[Any|str]): List of flags of the peak.
        analyte (Union[None, Analyte]): Analyte of the peak.
        mass_spectrum (np.ndarray): Mass spectrum of the peak.
    '''

    mass_spectrum: np.ndarray

    __slots__ = 'mass_spectrum'

    def __init__(self, rt: float, height: float, width: float, boarders: np.ndarray, mass_spectrum: np.ndarray,
                 flag:str|None=None, area=None):
        super().__init__(rt, height, width, boarders, area=area)
        self.mass_spectrum = mass_spectrum
        if flag:
            self.flags = [flag]

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

    def save_mass_spectrum(self, path:str) -> None:
        '''
        Saves the mass spectrum of the peak to a CSV file.
        Args:
            path (str): Path to the file where the mass spectrum will be saved.
        '''

        mass_spectrum = pd.DataFrame(self.mass_spectrum)
        mass_spectrum.to_csv(path, index=False)


