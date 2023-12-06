from pathlib import Path

import xarray as xr
from rdkit import Chem
from rdkit.Chem import Descriptors

from reaction.layout import Combinatorial_Layout
from reaction.transformation import Transformation


class Well_Plate(Combinatorial_Layout):

    '''
    A class for storing information about a combinatorial reaction layout in a well plate.

    Attributes:
        design (xr.DataArray): xarray DataArray containing the combinatorial design of the well plate.
        rows (dict): Dictionary containing the mapping of the well plate rows to the combinatorial dimensions of the
        well plate.
        columns (dict): Dictionary containing the mapping of the well plate columns to the combinatorial dimensions of
        the well plate.
    '''

    def __init__(self, layout_file:Path|str, transformation:Transformation, meta_data_file:str|None=None):
        super().__init__(layout_file, meta_data_file, transformation)
        self.design = xr.DataArray(self.array, dims=['x', 'y'], coords={'x': list(map(chr, range(65, 65+self.array.shape[0]))), 'y': list(range(1, self.array.shape[1]+1))}, attrs=self.meta_data)
        self.rows = {self.design.x.values[i]: self.layout.x.dropna().to_numpy(dtype=str)[i] for i in range(len(self.design.x.values))}
        self.columns = {self.design.y.values[i]: self.layout.y.dropna().to_numpy(dtype=str)[i] for i in range(len(self.design.y.values))}
        if meta_data_file:
            self.__extend_metadata()


    def __getitem__(self, pos:str) -> list|str:

        '''
        Takes a well position and returns the substrates in this position.
        '''

        item = self.design.loc[pos[0], int(pos[1:])].item()
        if '.' in item:
            item = item.split('.')
        return item

    def get_product_mw(self, pos:str):

        '''
        Returns the molecular weight of the product in the given position.

        Args:
            pos (str): Well position.

        Returns:
            float: Molecular weight of the product.

        '''

        product = self.get_product(pos)
        mol = Chem.MolFromSmiles(product)
        mw = Descriptors.ExactMolWt(mol)
        return mw

    def get_product_mz(self, pos:str):

        '''
        Returns the mass to charge ratio of the single charged product in the given position.
        Args:
            pos (str): Well position.

        Returns:
            float: Mass to charge ratio of the single charged product.
        '''

        mw = self.get_product_mw(pos)
        return round(mw, 0)


    def get_product(self, pos):

        '''
        Returns the product of the reaction in the given position.

        Args:
            pos (str): Well position.

        Returns:
            str: SMILES string of the product.

        '''

        subst = self[pos]
        product = self.transformation(subst)
        return product

    def get_substrate(self, pos, index=0):
        subst = self[pos]
        if isinstance(subst, list):
            subst = subst[index]
        return subst

    def __extend_metadata(self):

        '''
        Extends the metadata of the layout by replacing the '$row$' and '$column$' entries with the actual compounds
        from the layout.
        '''

        new_entries = []
        del_entries = []
        for i, stock in enumerate(self.meta_data['stock_solutions']):
            if stock['compound'] == '$row$':
                for well in stock['wells']:
                    new_entry = stock.copy()
                    new_entry['compound'] = self.rows[well]
                    new_entry['wells'] = [well]
                    new_entries.append(new_entry)
                    del_entries.append(i)
            if stock['compound'] == '$column$':
                for well in stock['wells']:
                    new_entry = stock.copy()
                    new_entry['compound'] = self.columns[int(well)]
                    new_entry['wells'] = [well]
                    new_entries.append(new_entry)
                    del_entries.append(i)
        for i in sorted(set(del_entries), reverse=True):
            del self.meta_data['stock_solutions'][i]
        self.meta_data['stock_solutions'].extend(new_entries)
