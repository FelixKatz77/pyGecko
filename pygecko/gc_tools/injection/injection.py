import _pickle as cPickle

import numpy as np

from pygecko.gc_tools.peak import Peak
from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.utilities import Utilities
from pygecko.visualization import Visualization




class Injection:

    '''
    Base class for all injections.

    Attributes:
        acq_method (str): Name of the acquisition method.
        instrument_name (str): Name of the instrument.
        sample_description (str): Description of the sample.
        sample_name (str): Name of the sample.
        sample_type (str): Type of the sample.
        vial_pos (int): Position of the sample's vial in the autosampler.
        internal_standard (Analyte): Internal standard of the sample.
        peaks (dict[float, Peak]): Peaks of the sample.
    '''

    acq_method: str
    instrument_name: str
    sample_description: str
    sample_name: str
    sample_type: str
    vial_pos: int
    internal_standard: Analyte|None
    peaks: dict[float, Peak]|None
    detector: None|str
    plate_pos: str | None
    analysis_settings: None
    chromatogram: None


    __slots__ = 'acq_method', 'instrument_name', 'sample_description', 'sample_name', 'sample_type', 'vial_pos', 'internal_standard', 'peaks', 'detector', 'plate_pos', 'analysis_settings', 'chromatogram'

    def __init__(self, metadata:dict, peaks:dict[float, Peak]|None=None, pos:bool=False):
        self.acq_method = metadata.get('AcqMethodName')
        self.instrument_name = metadata.get('InstrumentName')
        self.sample_description = metadata.get('SampleDescription')
        self.sample_name = metadata.get('SampleName')
        self.sample_type = metadata.get('SampleType')
        self.vial_pos = metadata.get('VialNumber')
        self.internal_standard = None
        self.peaks = peaks
        self.detector = None
        if pos:
            self.plate_pos = self.sample_name.split('-')[-1]
        else:
            self.plate_pos = None
        self.analysis_settings = None
        self.chromatogram = None

    def __getitem__(self, rt:float) -> Peak:

        '''
        Returns the peak with the given retention time.
        '''

        return self.peaks[rt]

    def __iter__(self) -> iter:

        '''
        Returns an iterator over the peaks.
        '''

        return iter(self.peaks.values())

    def __str__(self)-> str:

        '''
        Returns a string representation of the injection.
        '''

        if self.peaks is None:
            peak_count = 0
        else:
            peak_count = len(self.peaks)
        return f'{self.sample_name}: {peak_count} Peaks. {self.detector} Detection.'


    def set_internal_standard(self, rt:float|int, tolerance:float=0.05, name:str=None, smiles:str=None) -> None:

        '''
        Assigns the internal standard of the injection to the corresponding peak by creating an Analyte object for the
        internal standard and setting it as the peak's analyte.

        Args:
            rt (float|int): Retention time of the internal standard.
            tolerance (float): Tolerance for the retention time matching. Default is 0.05.
            name (str): Name of the internal standard. Default is None.
            smiles (str): SMILES string of the internal standard. Default is None.
        '''

        peak = self.flag_peak(rt, flag='standard', tolerance=tolerance)
        if not peak:
            raise ValueError(f'Error while setting internal standard: No peak found with retention time {rt} within tolerance {tolerance}.')
        self.internal_standard = Analyte(peak.rt, name=name, smiles=smiles)
        peak.analyte = self.internal_standard

    def flag_peak(self, rt: float, flag: str|None = None, tolerance: float = 0.05,
                  analyte: Analyte|None = None) -> None|Peak:

        '''
        Assigns a flag to the peak with the closest retention time to the given retention time within the tolerance.

        Args:
            rt (float): Retention time of the peak.
            flag (str|None): Flag to be assigned to the peak. Default is None.
            tolerance (float): Tolerance for the retention time matching. Default is 0.05.
            analyte (Analyte|None): Analyte object to be assigned to the peak. Default is None.

        Returns:
            None|Peak: The peak with the closest retention time to the given retention time or None if
            no peak was found within the tolerance.
        '''

        candidates = {}
        for peak in self.peaks.values():
            if Utilities.check_interval(rt, peak.rt, tolerance):
                deviation = (abs(rt - peak.rt))
                candidates[deviation] = peak
        if candidates:
            best_candidate = candidates[min(candidates)]
            best_candidate.flags.append(flag)
            best_candidate.analyte = analyte
            return best_candidate
        else:
            return None

    def match_ri(self, ri:float, tolerance:int=20, analyte:str|None=None, return_candidates:bool=False) -> Peak|None:

        '''
        Returns the peak with the closest retention index to the given retention index within the tolerance. Returns
        None if no peak was found within the tolerance.

        Args:
            ri (float): Retention index to match.
            tolerance(int): Tolerance for the retention index matching. Default is 20.
            analyte(str|None): Analyte object to be assigned to the peak. Default is None.

        Returns:
            Peak|None: The peak with the closest retention index to the given retention index or None if
             no peak was found within the tolerance.
        '''


        candidates = {}
        for peak in self.peaks.values():
            if peak.ri:
                if Utilities.check_interval(peak.ri, ri, tolerance):
                    deviation = abs(ri-peak.ri)
                    candidates[deviation] = peak
        if candidates:
            if return_candidates:
                return candidates
            else:
                if len(candidates) > 1:
                    for key in list(candidates.keys()):
                        if 'standard' in candidates[key].flags:
                            del candidates[key]
                peak = candidates[min(candidates)]
                if analyte:
                    peak.analyte = analyte
                return peak
        else:
            return None

    def get_plate_position(self):

        '''
        Returns the plate position of the sample as string.
        '''

        if self.plate_pos:
            return self.plate_pos
        else:
            raise AttributeError('Plate position is not assigned.')

    def set_plate_position(self, pos:str):

        '''
        Takes in a plate position and sets it as the sample's position.
        '''

        self.plate_pos = pos

    def view_chromatogram(self, path:str|None=None, **kwargs) -> None:

        '''
        Plots the chromatogram of the injection.

        Args:
            **kwargs: Keyword arguments for the visualization.
        '''

        Visualization.view_chromatogram(self, path=path, **kwargs)

    def _check_for_peak(self, chromatogram_slice) -> list[bool]:

        '''
        Takes in a chromatogram slice, returns a list of booleans indicating where a peak is present in the slice.
        '''

        bool_list = []

        for scan in chromatogram_slice:
            peak_assignment = False
            for peak in self.peaks.values():
                if peak.boarders[0] <= scan <= peak.boarders[1]:
                    peak_assignment = True
            bool_list.append(peak_assignment)
        return bool_list

    def save(self, filename: str) -> None:
        '''
        Saves an Injection to a .pkl file.

        Args:
            filename (str): Name of the file to save the injection to.
        '''

        with open(filename, 'wb') as outp:
            cPickle.dump(self, outp)

    def _check_for_missing_signal(self):

        '''Check if there is any part of the chromatogram with no signal and print the time ranges with no signal'''

        if np.any(self.chromatogram[1] == 0):
            print('No signal in the following time ranges:')
            ranges = Utilities.find_empty_ranges(self.chromatogram)
            for i, (start, end) in enumerate(Utilities.find_empty_ranges(self.chromatogram)):
                print(f'Range {i}: {start} - {end}')


def load_injection(filename) -> Injection:

    '''
    Loads an Injection from a .pkl file.

    Args:
        filename (str): Name of the file to load the ínjection from.

    Returns:
        Injection: Injection object loaded from the file.
    '''

    with open(filename, 'rb') as file:
        sequence = cPickle.load(file)
    return sequence





