from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from pygecko.gc_tools.peak.fid_peak import FID_Peak


class Quantification:

    '''
    A class wrapping functions for calculating yields based on peak areas.
    '''

    @staticmethod
    def quantify_ratio(peak:FID_Peak, standard:FID_Peak) -> int:

        '''
        Returns the relative quantity of an analyte based on the ratio of the analytes peak area to the peak area of a
        standard.

        Args:
            peak (FID_Peak): Peak to calculate the relative quantity for.
            standard (FID_Peak): Standard peak to use for quantification.

        Returns:
            yield_ (int): Relative quantity of analyte.
        '''

        yield_ = peak.area/standard.area
        return int(round(yield_, 0))

    @staticmethod
    def quantify_calibration(peak:FID_Peak, standard:FID_Peak, slope:float, intercept: float) -> int:

        '''
        Returns the relative quantity of an analyte based on a calibration curve.

        Args:
            peak (FID_Peak): Peak to calculate the relative quantity for.
            standard (FID_Peak): Standard peak to use for quantification.
            slope (float): Slope of the calibration curve.
            intercept (float): Intercept of the calibration curve.

        Returns:
            yield_ (int): Relative quantity of analyte.
        '''

        yield_ = slope * (peak.area/standard.area) + intercept
        return int(round(yield_, 0))

    @staticmethod
    def quantify_polyarc(peak:FID_Peak, std_peak:FID_Peak) -> int:

        '''
        Returns the relative quantity of an analyte based on the ratio of the analytes carbon-normalized peak area to
        the carbon-normalized peak area of a standard.

        Args:
            peak: Peak to calculate the relative quantity for.
            std_peak: Standard peak to use for quantification.

        Returns:
            yield_ (int): Relative quantity of analyte.
        '''

        analyte = peak.analyte.mol
        std = std_peak.analyte.mol
        if analyte:
            if std:
                C_analyte = Quantification.__get_c_count(analyte)
                C_std = Quantification.__get_c_count(std)
                yield_ = ((peak.area/C_analyte) / (std_peak.area/C_std)) * 100
                return int(round(yield_, 0))
            else:
                raise TypeError('No molecule assigned for Standard.')
        else:
            raise TypeError('No molecule assigned for Analyte.')


    @staticmethod
    def __get_c_count(mol:Mol):

        '''
        Takes in a molecule and returns the number of carbon atoms in the molecule.

        Args:
            mol (rdkit Mol object): Molecule to count the number of carbon atoms for.

        Returns:
            C_count (int): Number of carbon atoms in the molecule.
        '''

        C = Chem.MolFromSmarts('[#6]')
        C_count = len(mol.GetSubstructMatches(C))
        return C_count

