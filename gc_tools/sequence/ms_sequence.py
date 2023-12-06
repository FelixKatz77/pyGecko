from gc_tools.sequence.gc_sequence import GC_Sequence
from gc_tools.injection.ms_injection import MS_Injection


class MS_Sequence(GC_Sequence):

    '''
    Class to represent MS sequences.
    '''

    injections: dict[str:MS_Injection]

    def __init__(self, metadata: dict, injections: dict[str:MS_Injection]):
        super().__init__(metadata, injections)
        self.detector = 'MS'
