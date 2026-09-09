'''Fast regressions over small and already-tracked files from the pyGecko study.'''

from pathlib import Path

import numpy as np
import pytest

from pygecko.analysis import Analysis
from pygecko.gc_tools import FID_Sequence, MS_Sequence
from pygecko.parsers import FID_Base_Parser, MS_Base_Parser
from pygecko.reaction import Product_Array

from .conftest import FIXTURES


ROOT = Path(__file__).parents[2]
EXCERPTS = ROOT / 'tests/real_data/fixtures'
QUANTIFICATION = ROOT / 'notebooks/quantification'
PRODUCT = 'COC(=O)CCSc1c(Cl)ncn1C'


def test_xy_excerpt_detects_thiolation_product():
    injection = FID_Base_Parser.load_injection(
        EXCERPTS / 'thiolation_A1_product.xy', solvent_delay=0)

    injection.pick_peaks()

    peak = injection.peaks[6.920]
    assert peak.area > 0
    assert peak.boarders == pytest.approx([6.85865, 7.19931666666667])


def test_mzml_excerpt_detects_thiolation_product():
    injection = MS_Base_Parser.load_injection(
        EXCERPTS / 'thiolation_A1_product.mzML', pos=True)

    injection.pick_peaks()
    match = injection.match_mol(PRODUCT, check_iso=True)

    assert injection.sample_name == 'FKB-FA-060-A1'
    assert injection.scans.shape == (52, 174)
    assert match.rt == 5.941
    assert len(match.mass_spectrum) == 157


def test_thiolation_a1_reproduces_published_yield_and_retention_times(tmp_path):
    product_file = tmp_path / 'product.csv'
    product_file.write_text(f'product\n{PRODUCT}\n')
    layout = Product_Array(product_file)

    fid_injection = FID_Base_Parser.load_injection(
        QUANTIFICATION / 'FKB-FA-060-A1.xy', 2.7, pos=True)
    ms_injection = MS_Base_Parser.load_injection(
        FIXTURES / 'test_ri_calibration/FKB-FA-060-A1.mzML', pos=True)
    fid_injection.pick_peaks()
    ms_injection.pick_peaks()
    fid_injection.set_internal_standard(
        3.407, name='Dodecane', smiles='CCCCCCCCCCCC')
    ms_injection.set_internal_standard(
        2.260, name='Dodecane', smiles='CCCCCCCCCCCC')

    ms_calibration = MS_Base_Parser.load_ri_calibration(
        FIXTURES / 'test_ri_calibration/FKB-FA-060-RI.mzML', 12, rt=2.255)
    fid_calibration = FID_Base_Parser.load_ri_calibration(
        QUANTIFICATION / 'FKB-FA-060-RI.xy', 2.4, c_count=12, rt=3.394)
    ms_calibration.assign_ris(ms_injection)
    fid_calibration.assign_ris(fid_injection, alignment=True)

    ms_sequence = MS_Sequence({}, {ms_injection.sample_name: ms_injection})
    fid_sequence = FID_Sequence({}, {fid_injection.sample_name: fid_injection})
    result = Analysis.calc_plate_yield(ms_sequence, fid_sequence, layout)

    assert result.shape == (1, 1)
    assert result['quantity'][0, 0] == 75
    assert result['rt_ms'][0, 0] == 5.941
    assert result['rt_fid'][0, 0] == 6.920
    assert result['flags'][0, 0] == 0
