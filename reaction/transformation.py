from rdkit import Chem
from rdkit.Chem import AllChem

class Transformation:

    '''
    A class for storing information about a chemical transformation.

    Attributes:
        rxn_smarts (str): Reaction SMARTS string.
        transform (rdkit ChemicalReaction object): RDKit ChemicalReaction object.
    '''

    rxn_smarts: str
    transform: Chem.AllChem.ChemicalReaction

    __slots__ = 'rxn_smarts', 'transform'

    def __init__(self, rxn_smarts:str):
        self.rxn_smarts = rxn_smarts
        self.transform = AllChem.ReactionFromSmarts(rxn_smarts)

    def __call__(self, substrates) -> str:

        '''
        Returns the product of the transformation.

        Args:
            substrates (list[str]): List of SMILES strings of the substrates.

        Returns:
            product (str): SMILES string of the product.
        '''

        if len(substrates) != self.transform.GetNumReactantTemplates():
            print('Number of Substrates provided does not match transform!')
        else:
            mols = [Chem.MolFromSmiles(s) for s in substrates]
            #mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in substrates] #TODO: Fix this
            products = self.transform.RunReactants(mols)
            return Chem.MolToSmiles(products[0][0])

