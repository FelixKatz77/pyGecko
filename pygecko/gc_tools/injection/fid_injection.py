import numpy as np
import pandas as pd
from scipy import integrate
from scipy.signal import find_peaks

from pygecko.gc_tools.injection import Injection
from pygecko.gc_tools.peak import FID_Peak, Peak_Detection_FID
from pygecko.gc_tools.analysis import Analysis_Settings, Quantification
from pygecko.gc_tools.utilities import Utilities



class FID_Injection(Injection):

    '''
    Class to represent FID injections.

    Attributes:
        injector_pos (int): Position of the injector used for the injection.
        sample_number (int): Number of the sample in the sequence.
        acq_time (str): Acquisition time of the injection.
        data_method (Analysis_Settings): Data method used for the injection.
        solvent_delay (float): Solvent delay applied to the injection.
        chromatogram (np.ndarray): Chromatogram of the injection.
        processed_chromatogram (np.ndarray|None): Processed chromatogram of the injection.
        peaks (dict[float, FID_Peak]): Peaks of the injection.
        detector (str): Detector used for the injection.
    '''

    injector_pos: int
    sample_number: int
    acq_time: str
    analysis_settings: Analysis_Settings
    solvent_delay: float
    chromatogram: np.ndarray
    processed_chromatogram: np.ndarray|None
    peaks: dict[float, FID_Peak]|None
    detector: str

    __slots__ = 'injector_pos', 'sample_number', 'acq_time', 'analysis_settings', 'solvent_delay', 'chromatogram', 'processed_chromatogram', 'peaks', 'detector'

    peaks: None|list[FID_Peak]

    def __init__(self, metadata:dict, chromatogram:np.ndarray, solvent_delay:float|None=None, pos:bool=False):
        super().__init__(metadata, pos=pos)
        self.injector_pos = metadata.get('InjectorPosition')
        self.sample_number = metadata.get('SampleOrderNumber')
        self.acq_time = metadata.get('InjectionAcqDateTime')
        self.analysis_settings = Analysis_Settings(chromatogram)
        if solvent_delay: self.solvent_delay = solvent_delay
        else: self.solvent_delay = self.__set_solvent_delay(chromatogram)
        self.chromatogram = chromatogram[:,
                            Utilities.convert_time_to_scan(solvent_delay, self.analysis_settings.scan_rate):]
        self.processed_chromatogram = None
        self.peaks = None
        self.detector = 'FID'



    def baseline_correction(self, **kwargs:dict) -> None:

        '''
        Applies a baseline correction to the injection's chromatogram and sets the corrected chromatogram as the
        injection's processed chromatogram.

        Args:
            **kwargs: Keyword arguments for the baseline correction.
        '''

        self.analysis_settings.update(**kwargs)
        self.processed_chromatogram = Peak_Detection_FID.baseline_correction(self.chromatogram, self.analysis_settings)


    def pick_peaks(self, inplace: bool = True, **kwargs:dict) -> None|dict[float,FID_Peak]:

        '''
        Picks peaks from the injection's chromatogram.

        Args:
            inplace (bool): If True, the peaks are assigned to the injection's peaks attribute. Default is True.
            **kwargs: Keyword arguments for the peak picking.

        Returns:
            None|dict[float, FID_Peak]: The peaks of the injection if the inplace argument is False, None
            otherwise.
        '''

        self.analysis_settings.update(**kwargs)
        if not isinstance(self.processed_chromatogram, np.ndarray):
            self.baseline_correction()
        peaks = Peak_Detection_FID.pick_peaks(self.processed_chromatogram, self.analysis_settings)
        if inplace:
            self.peaks = peaks
        else:
            return peaks

    def integrate(self) -> None:

        '''
        Integrates the area under the curve of the injection's peaks and sets the area as the peak's area attribute.
        '''

        if self.peaks:
            for peak in self.peaks.values():
                area = integrate.simpson(self.chromatogram[1][round(peak.boarders[0]):round(peak.boarders[1])])
                peak.area = area
        else:
            print('Peaks list is empty.')

    def quantify(self, rt:float, method:str='polyarc', **kwargs) -> int:

        '''
        Returns the yield of the analyte with the given retention time calculated using the internal standard of the
        injection.

        Args:
            rt (float): Retention time of the analyte.
            method (str): Method to use for the quantification. Default is 'polyarc'.

        Returns:
            int: Yield of the analyte.
        '''

        if method == 'polyarc':
            yield_ = Quantification.quantify_polyarc(self.peaks[rt], self.peaks[self.internal_standard.rt])
            return yield_
        if method == 'calibration':
            yield_ = Quantification.quantify_calibration(self.peaks[rt], self.peaks[self.internal_standard.rt],
                                                         kwargs['slope'], kwargs['intercept'])
            return yield_

    def report(self, path:str) -> None:

        '''
        Writes a csv report for the injection to the given path.

        Args:
            path (str): Path to write the report to.
        '''

        injection_info = pd.DataFrame.from_dict({'Sample Name': [self.sample_name],
                                       'Sample Type': [self.sample_type],
                                       'Sample Description': [self.sample_description],
                                       'Acq. Method': [self.acq_method],
                                       'Injector Position': [self.injector_pos],
                                       'Acquisition Time': [self.acq_time.strftime('%d.%m.%Y; %H:%M:%S')],
                                       'Solvent Delay': [self.solvent_delay],
                                       'Detector': [self.detector],
                                       '': ''}, orient='index', columns=[0])
        peaks_info = pd.DataFrame.from_dict({i: ["{:.2f}".format(peak.rt), "{:.1f}".format(peak.area),
                                                 "{:.3f}".format(peak.width), "{:.2f}".format(peak.height)] for i, peak
                                             in enumerate(self.peaks.values(), start=1)},
                                            orient='index', columns=['RT [min]', 'Area', 'Width [min]', 'Height'])
        injection_info.to_csv(path, header=False)
        peaks_info.to_csv(path, mode='a')

    @staticmethod
    def __set_solvent_delay(chromatogram: np.ndarray) -> float:

        '''
        Returns the solvent delay of a chromatogram assuming the solvent peak is the highest peak.

        Args:
            chromatogram (np.ndarray): Chromatogram to detect the solvent delay for.

        Returns:
            float: Solvent delay for the chromatogram.
        '''

        peak_indices, peak_properties = find_peaks(chromatogram[1], width=0,
                                                   height=chromatogram[1].max(), rel_height=0.5)
        right_booarder = int(peak_properties['right_ips'][0].round(0))
        solvent_delay = chromatogram[0][right_booarder]
        return solvent_delay