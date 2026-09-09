'''Golden regressions over all three plates in Zenodo record 14316687.'''

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import pytest

from pygecko.analysis import Analysis
from pygecko.parsers import Agilent_FID_Parser, MS_Base_Parser
from pygecko.reaction import Reaction_Array, Transformation


pytestmark = [pytest.mark.realdata, pytest.mark.slow]
EXPECTED = Path(__file__).parent / 'expected'


@dataclass(frozen=True)
class PlateCase:
    name: str
    directory: str
    expected: str
    reaction: str
    fid_delay: float
    fid_standard: float
    ms_standard: float
    ms_calibration: str
    ms_calibration_carbon: int
    ms_calibration_rt: float
    fid_calibration_delay: float
    fid_calibration_rt: float
    ms_peak_kwargs: dict = field(default_factory=dict)


CASES = (
    PlateCase(
        'thiolation', 'Thiolation_Plate_GC_Data', 'thiolation_plate_yield.csv',
        '[c:1]1([a:2][a:3][a:4][a:5]1)[Cl,Br:6].[SH1:7][#6:8]>>'
        '[c:1]1([a:2][a:3][a:4][a:5]1)[S:7][#6:8]',
        2.7, 3.407, 2.260, 'FKB-FA-060-RI.mzML', 12, 2.255, 2.4, 3.394),
    PlateCase(
        'buchwald_hartwig', 'Buchwald-Hartwig_Plate_GC_Data', 'bh_plate_yield.csv',
        '[C,c:1][Nh1,Nh2,nh1:2].[Br,Cl:3][C,c:4]>>[C,c:1][N,n:2][C,c:4]',
        2.7, 4.593, 3.324, 'FBS-FA-033-RI-II.mzML', 10, 2.154, 2.7, 4.593,
        {'prominence_ms': 125}),
    PlateCase(
        'adhoc', 'AD-HoC_Plate_GC_Data', 'adhoc_plate_yield_paper.csv',
        '[C,c:1][Nh1,Nh2,nh1:2].[Br,Cl:3][C,c:4]>>[C,c:1][N,n:2][C,c:4]',
        3.0, 4.570, 3.320, 'FBS-FA-033-RI-II.mzML', 12, 3.320, 2.4, 4.570),
)


@pytest.fixture(scope='session')
def zenodo_data_root() -> Path:
    configured = os.environ.get('PYGECKO_REAL_DATA_DIR')
    if not configured:
        pytest.skip('set PYGECKO_REAL_DATA_DIR to extracted Zenodo record 14316687')
    root = Path(configured)
    if not root.is_dir():
        pytest.fail(f'PYGECKO_REAL_DATA_DIR is not a directory: {root}')
    return root


@pytest.mark.parametrize('case', CASES, ids=lambda case: case.name)
def test_plate_reproduces_published_results(case, zenodo_data_root, tmp_path):
    root = zenodo_data_root / case.directory
    layout = Reaction_Array(
        root / 'plate_layout.csv', Transformation(case.reaction),
        meta_data_file=root / 'meta_data.json')
    fid_sequence = Agilent_FID_Parser.load_sequence(root / 'FID', case.fid_delay, pos=True)
    ms_sequence = MS_Base_Parser.load_sequence(root / 'MS', pos=True)

    assert len(fid_sequence.injections) == 96
    assert len(ms_sequence.injections) == 96

    fid_sequence.pick_peaks()
    ms_sequence.pick_peaks(**case.ms_peak_kwargs)
    fid_sequence.set_internal_standard(
        case.fid_standard, name='Dodecane', smiles='CCCCCCCCCCCC')
    ms_sequence.set_internal_standard(
        case.ms_standard, name='Dodecane', smiles='CCCCCCCCCCCC')

    ms_calibration = MS_Base_Parser.load_ri_calibration(
        root / 'RI/MS' / case.ms_calibration,
        case.ms_calibration_carbon, rt=case.ms_calibration_rt)
    fid_calibration = Agilent_FID_Parser.load_ri_calibration(
        root / 'RI/FID', case.fid_calibration_delay,
        c_count=12, rt=case.fid_calibration_rt)
    ms_calibration.assign_ris(ms_sequence)
    fid_calibration.assign_ris(fid_sequence, alignment=True)

    output = tmp_path / f'{case.name}.csv'
    Analysis.calc_plate_yield(ms_sequence, fid_sequence, layout, path=output)
    actual = pd.read_csv(output, index_col=0)
    expected = pd.read_csv(EXPECTED / case.expected, index_col=0)

    if case.name == 'buchwald_hartwig':
        # The paper predates the overlap-border correction. Only this flagged peak changed.
        expected.loc['C9', 'Yield [%]'] = 67.0
        assert 'overlap' in fid_sequence.get_injection_by_pos('C9').peaks[7.997].flags

    pd.testing.assert_frame_equal(actual, expected, check_exact=True)
