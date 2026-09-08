import _pickle as cPickle
import json

from pygecko.gc_tools.history import Processing_Step, records_processing
from pygecko.gc_tools.peak import Peak
from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.utilities import Utilities




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
        history (list[Processing_Step]): Ordered record of the processing applied to the injection,
        starting with the parser's load call, so the state can be traced back to the raw data.
        _recording (bool): True while a recorded operation runs, so a nested call does not append a
        second step for work the caller never asked for.
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
    history: list[Processing_Step]
    _recording: bool


    __slots__ = 'acq_method', 'instrument_name', 'sample_description', 'sample_name', 'sample_type', 'vial_pos', 'internal_standard', 'peaks', 'detector', 'plate_pos', 'analysis_settings', 'chromatogram', 'history', '_recording'

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
        self.history = []
        self._recording = False

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


    @records_processing
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

    @records_processing
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

    @records_processing
    def match_ri(self, ri:float, tolerance:int=20, analyte:str|None=None,
                 return_candidates:bool=False) -> Peak|list[Peak]|None:

        '''
        Returns the peak with the closest retention index to the given retention index within the tolerance. Returns
        None if no peak was found within the tolerance.

        Args:
            ri (float): Retention index to match.
            tolerance(int): Tolerance for the retention index matching. Default is 20.
            analyte(str|None): Analyte object to be assigned to the peak. Default is None.
            return_candidates (bool): If True, returns every peak inside the tolerance window as a list
                instead of the single closest match, leaving the choice to the caller (Analysis
                disambiguates them by height ratio). Peaks flagged as the internal standard are kept in
                that list; they are only dropped from the single-match result. Default is False.

        Returns:
            Peak|list[Peak]|None: The peak with the closest retention index to the given retention index,
             the list of candidate peaks if return_candidates is True, or None if no peak was found within
             the tolerance.
        '''


        candidates = [peak for peak in self.peaks.values()
                      if peak.ri and Utilities.check_interval(peak.ri, ri, tolerance)]
        if candidates:
            if return_candidates:
                return candidates
            else:
                if len(candidates) > 1:
                    candidates = [peak for peak in candidates if 'standard' not in peak.flags]
                peak = min(candidates, key=lambda candidate: abs(ri - candidate.ri))
                if analyte:
                    peak.analyte = analyte
                return peak
        else:
            return None

    @records_processing
    def match_rt(self, rt:float, func=None, tolerance:float=1/60, analyte:str|None=None,
                 return_candidates:bool=False, exclude_standard:bool=True) -> Peak|list[Peak]|None:

        '''
        Returns the peak with the closest retention time to the (optionally transformed) target retention time
        within the tolerance. Returns None if no peak was found within the tolerance.

        This is the retention-time analogue of ``match_ri`` for split-GC data, where the FID and MS traces
        originate from a single injection and therefore share a retention-time axis (up to a small, near-constant
        splitter dead-volume offset). Instead of converting retention times to retention indices, the expected
        retention time in this injection's trace is obtained by applying ``func`` to the source (e.g. MS)
        retention time. ``func`` defaults to the identity ``y = x`` (assume identical retention times); a drift
        model (e.g. ``lambda x: a*x + b``) can be supplied later without changing the matching logic.

        Args:
            rt (float): Source retention time to match (e.g. the MS peak's retention time), in minutes.
            func (callable|None): Mapping from the source retention time to the expected retention time in this
                injection's trace. Defaults to the identity function ``lambda x: x``.
            tolerance (float): Half-width of the retention-time matching window, in minutes. Defaults to 1/60
                (one second).
            analyte (str|None): Analyte object to assign to the matched peak. Default is None.
            return_candidates (bool): If True, returns every peak inside the tolerance window as a list
                instead of the single closest match, leaving the choice to the caller (Analysis
                disambiguates them by height ratio). Default is False.
            exclude_standard (bool): If True, peaks flagged as the internal standard are never returned as a
                match. Default is True.

        Returns:
            Peak|list[Peak]|None: The peak with the closest retention time to the target, the list of
            candidate peaks if return_candidates is True, or None if no peak was found within the tolerance.
        '''

        if func is None:
            func = lambda x: x
        target = func(rt)
        candidates = [peak for peak in self.peaks.values()
                      if not (exclude_standard and 'standard' in peak.flags)
                      and Utilities.check_interval(peak.rt, target, tolerance)]
        if not candidates:
            return None
        if return_candidates:
            return candidates
        peak = min(candidates, key=lambda candidate: abs(candidate.rt - target))
        if analyte:
            peak.analyte = analyte
        return peak

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

        # Imported here, not at module scope: visualization imports gc_tools, so a
        # module-level import makes pygecko.visualization unimportable on its own.
        from pygecko.visualization import Visualization

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

    def record_step(self, operation:str, parameters:dict) -> None:

        '''
        Appends a processing step to the injection's history.

        For provenance no decorated injection method produces: the parser's load call, which is the
        first step and the point a replay would start from, and RI_Calibration.assign_ris, which
        mutates the injection's peaks from outside the injection. The operation name and its
        arguments arrive as plain data, so gc_tools still knows nothing about file formats or the
        calibration's internals.

        Args:
            operation (str): Qualified name of the operation as Class.method.
            parameters (dict): Arguments that reproduce the operation.
        '''

        self.history.append(Processing_Step(operation, parameters))

    def history_to_json(self, path:str|None=None) -> str:

        '''
        Returns the injection's processing history as a JSON document, writing it to path when one
        is given.

        Args:
            path (str|None): Path to write the document to. Default is None, which only returns it.

        Returns:
            str: The history as a JSON array of steps, each with its operation, parameters, resolved
            settings, timestamp and version.
        '''

        document = json.dumps([step.to_dict() for step in self.history], indent=2)
        if path:
            with open(path, 'w') as outp:
                outp.write(document)
        return document

    def __setstate__(self, state:tuple) -> None:

        '''
        Restores an Injection from its pickled state, defaulting attributes added since the file was
        written.

        Pickle files are coupled to the class layout, so injections saved before the processing
        history existed carry no history or _recording entry. Both are defaulted rather than
        reconstructed: the steps that produced those objects were never recorded, and an empty
        history is the honest statement of that.
        '''

        _, slots = state
        for name, value in (slots or {}).items():
            setattr(self, name, value)
        if not hasattr(self, 'history'):
            self.history = []
        if not hasattr(self, '_recording'):
            self._recording = False

    def save(self, filename: str) -> None:
        '''
        Saves an Injection to a .pkl file.

        Args:
            filename (str): Name of the file to save the injection to.
        '''

        with open(filename, 'wb') as outp:
            cPickle.dump(self, outp)

    def _check_for_missing_signal(self):

        '''Check the chromatogram for stretches with no signal and print a one-line summary.

        A real MS TIC carries dozens of short dropouts per injection, so the ranges are summarised
        rather than listed: itemising them buries every other message a plate run prints. Use
        Utilities.find_empty_ranges directly to inspect the individual gaps.
        '''

        ranges = Utilities.find_empty_ranges(self.chromatogram)
        if not ranges:
            return
        dead_time = sum(end - start for start, end in ranges)
        run_time = self.chromatogram[0][-1] - self.chromatogram[0][0]
        start, end = max(ranges, key=lambda r: r[1] - r[0])
        print(f'{self.sample_name}: {len(ranges)} empty ranges, {dead_time:.2f} min dead of '
              f'{run_time:.2f} min, longest {end - start:.3f} min at {start:.2f}-{end:.2f}')


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





