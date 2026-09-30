from typing import Literal

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from pygecko.gc_tools.chromatogram import Chromatogram
from pygecko.gc_tools.history import records_processing
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
        chromatogram (Chromatogram): FID trace of the injection from the solvent delay on; its processed
        signal is filled by the baseline correction.
        peaks (dict[float, FID_Peak]): Peaks of the injection.
        detector (str): Detector used for the injection.
        _baseline_settings (tuple|None): The baseline settings chromatogram.processed was computed with, or
        None if it has not been computed.
    '''

    injector_pos: int
    sample_number: int
    acq_time: str
    analysis_settings: Analysis_Settings
    solvent_delay: float
    chromatogram: Chromatogram
    peaks: dict[float, FID_Peak]|None
    detector: str
    _baseline_settings: tuple|None

    __slots__ = 'injector_pos', 'sample_number', 'acq_time', 'analysis_settings', 'solvent_delay', 'chromatogram', 'peaks', 'detector', '_baseline_settings'

    peaks: None|list[FID_Peak]

    def __init__(self, metadata:dict, chromatogram:Chromatogram, solvent_delay:float|None=None, pos:bool=False):
        super().__init__(metadata, pos=pos)
        self.injector_pos = metadata.get('InjectorPosition')
        self.sample_number = metadata.get('SampleOrderNumber')
        self.acq_time = metadata.get('InjectionAcqDateTime')
        self.analysis_settings = Analysis_Settings()
        if solvent_delay is not None:
            self.solvent_delay = solvent_delay
        else:
            self.solvent_delay = self.__set_solvent_delay(chromatogram)
        self.chromatogram = chromatogram[Utilities.convert_time_to_scan(self.solvent_delay, chromatogram.scan_rate):]
        self.peaks = None
        self.detector = 'FID'
        self._baseline_settings = None

    def __setstate__(self, state:tuple) -> None:

        '''
        Restores an FID_Injection from its pickled state, defaulting the baseline-settings cache key
        that files written before it existed do not carry, so the next pick recomputes the baseline.
        '''

        super().__setstate__(state)
        if not hasattr(self, '_baseline_settings'):
            self._baseline_settings = None



    @records_processing
    def baseline_correction(self, **kwargs:dict) -> None:

        '''
        Applies a baseline correction to the injection's chromatogram and sets the corrected chromatogram as the
        injection's processed chromatogram.

        Args:
            **kwargs: Keyword arguments for the baseline correction.
        '''

        self.analysis_settings.update(**kwargs)
        self.chromatogram.processed = Peak_Detection_FID.baseline_correction(self.chromatogram, self.analysis_settings)
        self._baseline_settings = self.__baseline_settings()


    @records_processing
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
        # The processed signal depends on the baseline settings only, so it is recomputed when they
        # change and reused when only the detection settings, time_range among them, do.
        if self._baseline_settings != self.__baseline_settings():
            self.baseline_correction()
        peaks = Peak_Detection_FID.pick_peaks(self.chromatogram, self.analysis_settings)
        if inplace:
            self.peaks = peaks
        else:
            return peaks

    @records_processing
    def integrate(self) -> None:

        '''
        Integrates the area under the curve of the injection's peaks between their boarders and sets
        the area as the peak's area attribute.

        Integrates the baseline corrected chromatogram, applying the baseline correction first if it
        has not been applied yet, so the areas agree with the ones pick_peaks already computed.

        Raises:
            ValueError: If the injection has no peaks.
        '''

        if not self.peaks:
            raise ValueError(f'{self.sample_name}: no peaks to integrate; call pick_peaks first.')
        if self.chromatogram.processed is None:
            self.baseline_correction()
        for peak in self.peaks.values():
            # Boarders are retention times in minutes, not scan indices, so they are looked up
            # on the chromatogram's own time axis - exactly, because pick_peaks took them from
            # that same axis. The signal is the baseline-corrected one pick_peaks integrated:
            # quantification divides one area by another and a baseline offset does not cancel
            # between peaks of different width.
            start, end = np.searchsorted(self.chromatogram.time, peak.boarders)
            peak.area = Peak_Detection_FID.peak_area(self.chromatogram.time, self.chromatogram.processed, start, end)

    @records_processing
    def quantify(self, rt:float, method:Literal['polyarc', 'ratio', 'calibration']='polyarc', **kwargs) -> float:

        '''
        Returns the yield of the analyte with the given retention time calculated using the internal standard of the
        injection.

        Args:
            rt (float): Retention time of the analyte, as the key it has in peaks.
            method (str): Method to use for the quantification: 'polyarc', 'ratio' or 'calibration'. Default is
                'polyarc'. 'calibration' takes slope and intercept as keyword arguments.

        Returns:
            float: Yield of the analyte in percent.

        Raises:
            ValueError: If no internal standard is set or the method is unknown.
            KeyError: If no peak has the retention time rt.
        '''

        if self.internal_standard is None:
            raise ValueError(f'{self.sample_name}: no internal standard set.')
        if rt not in self.peaks:
            raise KeyError(f'{self.sample_name}: no peak at {rt} min; peaks are keyed by retention time '
                           f'rounded to 3 decimals.')
        peak, standard = self.peaks[rt], self.peaks[self.internal_standard.rt]
        if method == 'polyarc':
            return Quantification.quantify_polyarc(peak, standard)
        if method == 'ratio':
            return Quantification.quantify_ratio(peak, standard)
        if method == 'calibration':
            return Quantification.quantify_calibration(peak, standard, kwargs['slope'], kwargs['intercept'])
        raise ValueError(f'Unknown quantification method {method!r}; expected polyarc, ratio or calibration.')

    def __baseline_settings(self) -> tuple:

        '''
        Returns the configured settings chromatogram.processed depends on; time_range is not one of them.
        '''

        return self.analysis_settings.savgol_window, self.analysis_settings.max_half_window

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
    def __set_solvent_delay(chromatogram: Chromatogram) -> float:

        '''
        Returns the solvent delay of a chromatogram assuming the solvent peak is the highest peak.

        Args:
            chromatogram (Chromatogram): Chromatogram to detect the solvent delay for.

        Returns:
            float: Solvent delay for the chromatogram.
        '''

        peak_indices, peak_properties = find_peaks(chromatogram.intensity, width=0,
                                                   height=chromatogram.intensity.max(), rel_height=0.5)
        right_booarder = int(peak_properties['right_ips'][0].round(0))
        solvent_delay = chromatogram.time[right_booarder]
        return solvent_delay

