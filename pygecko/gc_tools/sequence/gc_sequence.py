import _pickle as cPickle
import pickle as pickle
from pygecko.gc_tools.injection.injection import Injection
from pygecko.gc_tools.analyte import Analyte
from pathlib import Path


class GC_Sequence:

    '''
    Base class for all GC sequences.

    Attributes:
        sequence_name (str): Name of the sequence.
        instrument_name (str): Name of the instrument the sequence was measured on.
        injections (dict[str, Injection]): Injections of the sequence.
        internal_standard (Analyte|None): Internal standard of the sequence.
    '''

    sequence_name: str
    instrument_name: str
    injections: dict[str, Injection]
    internal_standard: Analyte|None
    detector: None|str

    __slots__ = 'sequence_name', 'instrument_name', 'injections', 'internal_standard', 'detector'

    def __init__(self, metadata: dict, injections: dict[str:Injection]):
        self.sequence_name = metadata.get('sequence_name')
        self.instrument_name = metadata.get('instrument_name')
        self.injections = injections
        self.internal_standard = None
        self.detector = None

    def __getitem__(self, sample_name:str) -> Injection:

        '''
        Returns the injection with the given sample name.
        '''

        return self.injections[sample_name]

    def __iter__(self) -> iter:

        '''
        Returns an iterator over the injections.
        '''

        return iter(self.injections.values())

    def __len__(self)-> int:

        '''
        Returns the number of injections.
        '''

        return len(self.injections)

    def __contains__(self, sample_name:str) -> bool:

        '''
        Returns True if the sequence contains an injection with the given sample name and False otherwise.
        '''

        return sample_name in self.injections

    def __str__(self) -> str:

        '''
        Returns a string representation of the sequence.
        '''

        return f'{self.sequence_name}: {len(self)} Injections, {self.detector} Detection.'

    def __delete__(self, injection_name:str)-> None:

        '''
        Deletes the injection with the given name from the sequence.
        '''

        del self.injections[injection_name]

    def pick_peaks(self, **kwargs):

        '''
        Picks peaks from the injections' chromatograms.
        Args:
            **kwargs: Keyword arguments for the peak picking.
        '''

        for injection in self.injections.values():
            injection.pick_peaks(**kwargs)

    def set_internal_standard(self, rt:float|int, tolerance:float=0.05, name:str=None, smiles:str=None) -> None:

        '''
        Assigns the internal standard of the sequence to the corresponding peaks by creating an Analyte object for the
        internal standard and setting it as the peaks' analyte.

        Args:
            rt (float|int): Retention time of the internal standard.
            tolerance (float): Tolerance for the retention time matching. Default is 0.05.
            name (str): Name of the internal standard. Default is None.
            smiles (str): SMILES string of the internal standard. Default is None.
        '''

        self.internal_standard = Analyte(rt, name=name, smiles=smiles)
        for injection in self.injections.values():
            injection.set_internal_standard(rt, tolerance=tolerance, name=name, smiles=smiles)

    def save(self, filename: str) -> None:

        '''
        Saves a GC sequence to a .pkl file.

        Args:
            filename (str): Name of the file to save the sequence to.
        '''

        with open(filename, 'wb') as outp:
            cPickle.dump(self, outp)

def load_sequence(filename:str) -> GC_Sequence:

    '''
    Loads a GC sequence from a .pkl file.

    Args:
        filename (str): Name of the file to load the sequence from.

    Returns:
        GC_Sequence: GC_Sequence object loaded from the file.
    '''
    load_path = Path(filename)
    with open(load_path, 'rb') as file:
        sequence = pickle.load(file)
    return sequence

def save_sequence(sequence:GC_Sequence, filename:str) -> None:
    '''
    Saves a GC sequence to a .pkl file.

    Args:
        filename (str): Name of the file to save the sequence to.
    '''
    save_path = Path.cwd().joinpath(filename)
    with open(save_path, 'wb') as outp:
        pickle.dump(sequence, outp)


