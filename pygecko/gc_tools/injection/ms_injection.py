from collections import defaultdict

import numpy as np
import pandas as pd
from brainpy import isotopic_variants
from rdkit import Chem
from rdkit.Chem import Descriptors

from pygecko.gc_tools.analyte import Analyte
from pygecko.gc_tools.analysis import Analysis_Settings
from pygecko.gc_tools.injection import Injection
from pygecko.gc_tools.peak import Peak_Detection_MS, MS_Peak



class MS_Injection(Injection):

    '''
    Class to represent MS injections.

    Attributes:
        chromatogram (np.ndarray): Chromatogram of the injection.
        peaks (dict[float, MS_Peak]): Peaks of the injection.
        scans (pd.DataFrame): Scans of the injection.
        detector (str): Detector used for the injection.
        analysis_settings (Analysis_Settings): Data method used for the injection.
        solvent_delay (float): Solvent delay applied to the injection.
    '''

    chromatogram: np.ndarray
    peaks: dict[float, MS_Peak]|None
    scans: pd.DataFrame
    detector: str
    analysis_settings: Analysis_Settings
    solvent_delay: float

    __slots__ = 'chromatogram', 'peaks', 'scans', 'detector', 'analysis_settings', 'solvent_delay'

    def __init__(self, metadata:dict|None, chromatogram:np.ndarray, peaks:dict|None, scans:pd.DataFrame, pos:bool=False):
        super().__init__(metadata, pos=pos)
        self.chromatogram = chromatogram
        self.peaks = peaks
        self.scans = scans
        self.detector = 'MS'
        self.analysis_settings = Analysis_Settings(chromatogram)
        self.solvent_delay = chromatogram[0][0]

    def match_mz(self, mz:float) -> MS_Peak|list[MS_Peak]|None:

        '''
        Returns the peak or a list of peaks whose mass spectra contain the given m/z. Returns None if no peak was found
        matching the criteria.

        Args:
            mz (float): m/z to match.

        Returns:
            MS_Peak|list[MS_Peak]|None: The peak or a list of peaks whose mass spectra contain the given m/z or None if
            no peak was found matching the criteria.

        '''

        candidates = []
        for rt, peak in self.peaks.items():
            if mz in peak.mass_spectrum['mz']:
                index = np.where(peak.mass_spectrum['mz'] == mz)[0]
                if peak.mass_spectrum['rel_intensity'][index][0] > 4 and mz > peak.mass_spectrum['mz'].max()*(2/3):
                    candidates.append(peak)
        if candidates:
            if len(candidates) > 1:
                print(f'Multiple peaks with m/z {mz} fitting the calculated isotope pattern were found for {self.sample_name}.')
                return candidates
            else:
                return candidates[0]
        return None

    def match_mol(self, smiles:str) -> MS_Peak|None:

        '''
        Returns the peak with the lowest isotope error for the m/z corresponding to the given molecule's parent
        peak. Returns None if no peak was found matching the criteria.

        Args:
            smiles (str): SMILES string of the analyte.

        Returns:
            MS_Peak|None: The peak with the lowest isotope error for the given molecule or None if no peak was
            found matching the criteria.
        '''

        mol = Chem.MolFromSmiles(smiles)
        mz = round(Descriptors.ExactMolWt(mol), 0)
        peak = self.__match_mz_mol(mz, smiles=smiles)
        if peak:
            analyte = Analyte(peak.rt, smiles=smiles)
            peak.analyte = analyte
        return peak

    def __match_mz_mol(self, mz:float, smiles:str) -> MS_Peak|None:

        '''
        Returns the peak with the lowest isotope error for a given m/z and molecule. Returns None if no peak was found
        matching the criteria.

        Args:
            smiles (str): SMILES string of the analyte.
            mz (float): m/z to match.

        Returns:
            MS_Peak|None: The peak with the lowest isotope error for the given m/z and molecule or None if no peak was
            found matching the criteria.
        '''

        candidates = {}
        for rt, peak in self.peaks.items():
            if mz in peak.mass_spectrum['mz']:
                index = np.where(peak.mass_spectrum['mz'] == mz)[0]
                if peak.mass_spectrum['rel_intensity'][index][0] > 4 and mz > peak.mass_spectrum['mz'].max() * (
                        2 / 3):  # TODO: Check if this is a good decision.
                    isotope_error = self.__isotope_check(smiles, peak, mz)
                    if isotope_error:
                        candidates[isotope_error] = peak
        if candidates:
            if len(candidates) > 1:
                print(f'Multiple peaks with m/z {mz} fitting the calculated isotope pattern were found for {self.sample_name}.')
            peak = candidates[min(candidates)]
            peak.analyte = Analyte(peak.rt, smiles=smiles)
            return peak
        return None
    def pick_peaks(self, inplace: bool = True, **kwargs: dict) -> None|dict[float, MS_Peak]:

        '''
        Picks peaks from the injection's chromatogram.

        Args:
            inplace (bool): If True, the peaks are assigned to the injection's peaks attribute. Default is True.
            **kwargs: Keyword arguments for the peak picking.

        Returns:
            None|dict[float, MS_Peak]: The peaks of the injection if the inplace argument is False, None
            otherwise.
        '''

        self.analysis_settings.update(**kwargs)
        peaks = Peak_Detection_MS.pick_peaks(self.chromatogram, self.scans, self.analysis_settings)
        if inplace:
            self.peaks = peaks
        else:
            return peaks

    def __isotope_check(self, smiles:str, peak:MS_Peak, mz:float) -> float|None:

        '''
        Returns True if the isotope peak of a given m/z is present in the peak's mass spectrum and False otherwise.

        Args:
            smiles (str): SMILES string of the analyte.
            peak (MS_Peak): Peak to check.
            mz (float): m/z to check.

        Returns:
            float|None: The difference between the theoretical and the measured ratio of the isotope peak of a given m/z
            if the difference is smaller than 5% and None otherwise.
        '''

        if 'Cl' in smiles or 'Br' in smiles:
            diff = 2
        else:
            diff = 1
        return self.__isotopic_ratio_check(smiles, peak, mz, diff)


    def __isotopic_ratio_check(self,smiles:str, peak:MS_Peak, mz:float, diff:int) -> float|None:

        '''
        Returns True if the ratio of the isotope peak of a given m/z is within 6% of the theoretical ratio and False if
        otherwise.

        Args:
            smiles (str): SMILES string of the analyte.
            peak (MS_Peak): Peak to check.
            mz (float): m/z to check.
            diff (int): Difference between the m/z of the isotope peak and the m/z of the peak to check.

        Returns:
            float|None: The difference between the theoretical and the measured ratio of the isotope peak of a given m/z
            if the difference is smaller than 5% and None otherwise.
        '''

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        i = np.where(peak.mass_spectrum['mz'] == mz)[0]
        j = np.where(peak.mass_spectrum['mz'] == mz+diff)[0]
        if not i or not j:
            return False
        ratio = peak.mass_spectrum['rel_intensity'][j][0]/peak.mass_spectrum['rel_intensity'][i][0]
        mol_formula = self.__get_mol_formula_dict(mol)
        isotopic_dist = isotopic_variants(mol_formula, npeaks=3, charge=0)
        theo_ratio = isotopic_dist[diff].intensity/isotopic_dist[0].intensity
        if (theo_ratio - 0.055) < ratio < (theo_ratio + 0.055):
            return abs(theo_ratio - ratio)
        else:
            return None

    def __get_mol_formula_dict(self, mol) -> dict:

        '''
        Returns a dictionary with the molecule's elements as keys and the number of atoms of each element as values.

        Args:
            mol (rdkit Mol Object): Molecule to get the formula of.

        Returns:
            dict: Dictionary with the molecule's elements as keys and the number of atoms of each element as values.

        '''

        mol_formula = defaultdict(lambda : 0)
        for atom in mol.GetAtoms():
            mol_formula[atom.GetSymbol()] += 1
        return mol_formula
