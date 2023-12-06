from rdkit import Chem
from rdkit.Chem import Descriptors

class Analyte:

    '''
    Analyte class for storing information about analytes.

    Attributes:
        rt (float): Retention time of analyte.
        name (str): Name of analyte.
        smiles (str): SMILES string of analyte.
        mol (Mol): RDKit molecule object of analyte.
    '''

    rt: float|int
    name: str
    smiles: str
    mz: None|int

    __slots__ = 'rt', 'name', 'smiles', 'mol', 'mz'

    def __init__(self, rt:float|int, name:str=None, smiles:str=None):
        self.rt = rt
        self.name = name
        self.smiles = smiles
        if smiles:
            self.mol = Chem.MolFromSmiles(smiles)
            self.mz = round(Descriptors.ExactMolWt(self.mol), 0)
        else:
            self.mol = None
            self.mz = None

    def __str__(self) -> str:
        '''
        Returns a string representation of the analyte.
        '''
        return f'{self.name} ({self.smiles})'



