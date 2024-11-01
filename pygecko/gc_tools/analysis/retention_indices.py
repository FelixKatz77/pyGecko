import warnings
from bisect import bisect_left
import numpy as np
from scipy import stats
from pygecko.gc_tools.injection import Injection, FID_Injection, MS_Injection
from pygecko.gc_tools.sequence import GC_Sequence
from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.peak import Peak


class RI_Calibration:

    '''
    A class for storing information about a retention index calibration and calculating retention indices.

    Attributes:
        calibration (FID_Injection|MS_Injection): Calibration injection.
        alkanes (np.ndarray): Array of alkanes with columns for smiles, carbon count and retention time.
        gradient (float): Gradient of the linear fit of the calibration.
        intercept (float): Intercept of the linear fit of the calibration.
    '''

    calibration: FID_Injection|MS_Injection
    alkanes: np.ndarray
    gradient: float
    intercept: float

    __slots__ = 'calibration', 'alkanes', 'gradient', 'intercept'

    def __init__(self, injection:FID_Injection|MS_Injection, c_count:int, rt:float):
        self.calibration = injection
        self.calibration.pick_peaks()
        self.__identify_alkanes(c_count, rt)
        self.alkanes = self.__construct_alkanes_array()
        self.gradient, self.intercept = self.__fit_ri_calibration()

    def calculate_ri_fit(self, rt):

        '''
        Returns the retention index for a given retention time calculated using the linear fit of the calibration.

        Args:
            rt (float): Retention time.

        Returns:
            ri (float): Retention index.
        '''

        ri = self.gradient * rt + self.intercept
        return ri

    def calculate_ri(self, peak:Peak, align_factor=0) -> float|None:

        '''
        Returns the retention index for a given peak.

        Args:
            peak (Peak): Peak to calculate the retention index for.
            align_factor (float): Alignment factor for the retention time. Default is 0.

        Returns:
            ri (float): Retention index.
        '''

        calibration_rts = self.alkanes['rt']
        calibration_c = self.alkanes['c_count']
        i = bisect_left(calibration_rts, peak.rt)
        if 0 < i < (len(calibration_rts)-1):
            alkane_rts = (calibration_rts[i-1], calibration_rts[i])
            c_counts = (calibration_c[i-1], calibration_c[i])
            ri = 100 * (c_counts[0] + ((peak.rt - (alkane_rts[0])) / (
                    (alkane_rts[1]) - (alkane_rts[0]))))
        else:
            gradient, intercept = self.__fit_ri_calibration(align_factor)
            ri = gradient * peak.rt + intercept
        return ri

    def assign_ris(self, gc_data: Injection|GC_Sequence, alignment=False) -> None:

        '''
        Assigns retention indices to peaks in a GC_Sequence or Injection.

        Args:
            gc_data(Injection|GC_Sequence): GC_Sequence or Injection to assign retention indices to.
            alignment (bool): If True, retention times are aligned to the internal standard. Default is False.
        '''

        if isinstance(gc_data, Injection):
            self.__assign_ris_injection(gc_data, alignment=alignment)
        elif isinstance(gc_data, GC_Sequence):
            self.__assign_ris_sequence(gc_data, alignment=alignment)
        else:
            raise TypeError(f'{type(gc_data)} is not supported by {self.assign_ris}.')

    def __assign_ris_injection(self, injection:Injection, alignment=False):

        '''
        Assigns retention indices to peaks in an Injection.

        Args:
            injection (Injection): Injection to assign retention indices to.
            alignment (bool): If True, retention times are aligned to the internal standard. Default is False.
        '''

        for rt, peak in injection.peaks.items():
            if alignment:
                align_factor = self.__calculate_align_factor(injection.internal_standard)
            else:
                align_factor = 0
            ri = self.calculate_ri(peak, align_factor=align_factor)
            if ri:
                peak.ri = ri
            else:
                warnings.warn(f'Could not assign RI for {injection.detector}-Peak at {rt} min in {injection.sample_name}. RT out of calibration range.')

    def __assign_ris_sequence(self, sequence:GC_Sequence, alignment=False):

        '''
        Assigns retention indices to peaks in a GC_Sequence.

        Args:
            sequence (GC_Sequence): GC_Sequence to assign retention indices to.
            alignment (bool): If True, retention times are aligned to the internal standard. Default is False.
        '''

        for injection in sequence.injections.values():
            self.__assign_ris_injection(injection, alignment=alignment)

    def __identify_alkanes(self, c_count:int, rt:float):

        '''
        Identifies alkanes in the calibration injection based on the carbon count and retention time given.

        Args:
            c_count (int): Carbon count of an alkane present in the calibration.
            rt (float): Retention time of the corresponding alkane.
        '''

        seed_peak = self.calibration.flag_peak(rt, flag=f'C{c_count}')
        seed_rt = seed_peak.rt
        rts = list(self.calibration.peaks.keys())
        index = rts.index(seed_rt)
        for i in range(index, len(self.calibration.peaks)):
            self.calibration.peaks[rts[i]].analyte = Analyte(rts[i], smiles='C'*(c_count + i - index))
            self.calibration.peaks[rts[i]].flags.append('Alkane_RI')
        for i in range(0, index):
            self.calibration.peaks[rts[i]].analyte = Analyte(rts[i], smiles='C'*(c_count + i - index))
            self.calibration.peaks[rts[i]].flags.append('Alkane_RI')

    def __construct_alkanes_array(self) -> dict[float:int]:

        '''
        Returns an array of alkanes with columns for smiles, carbon count and retention time.

        Returns:
            alkanes (np.ndarray): Array of alkanes with columns for smiles, carbon count and retention time.
        '''

        _alkanes_list = []
        for rt, peak in self.calibration.peaks.items():
            if 'Alkane_RI' in peak.flags:
                _alkanes_list.append((peak.analyte.smiles, len(peak.analyte.smiles), rt))
            else:
                pass
        alkanes = np.array(_alkanes_list, dtype=[('smiles', 'U50'), ('c_count', int), ('rt', float)])
        return alkanes

    def __fit_ri_calibration(self, align_factor:float=0):

        '''
        Returns the gradient and intercept of the linear fit of the calibration.
        Args:
            align_factor (float): Alignment factor for the retention time. Default is 0.

        Returns:
            gradient (float): Gradient of the linear fit of the calibration.
            intercept (float): Intercept of the linear fit of the calibration.
        '''

        rts = []
        ris = []
        for i, c_count in enumerate(self.alkanes['c_count']):
            rt = self.alkanes['rt'][i]
            ris.append(100 * c_count)
            rts.append(rt+align_factor)
        gradient, intercept, r_value, p_value, slope_std_error = stats.linregress(rts, ris)
        return gradient, intercept

    def __calculate_align_factor(self, internal_standard:Analyte):

        '''
        Returns the alignment factor for the retention time of the internal standard.

        Args:
            internal_standard (Analyte): Internal standard to align the retention time to.

        Returns:
            align_factor (float): Alignment factor for the retention time.
        '''

        index = np.where(self.alkanes['smiles'] == internal_standard.smiles)[0]
        align_factor = internal_standard.rt - self.alkanes['rt'][index][0]
        return align_factor






