from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.sequence.gc_sequence import GC_Sequence



class FID_Sequence(GC_Sequence):

    '''
    Class to represent FID sequences.
    '''

    injections: dict[str:FID_Injection]

    def __init__(self, metadata:dict|None, injections:dict[str, FID_Injection]):
        super().__init__(metadata, injections)
        self.detector = 'FID'
