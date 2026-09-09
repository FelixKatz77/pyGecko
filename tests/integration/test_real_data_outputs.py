'''Output-format regressions driven by the attributed thiolation study metadata.'''

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pygecko.data_handling import PDF_Report
from pygecko.data_handling import reports
from pygecko.reaction import Reaction_Array, Reaction_Parser, Transformation


FIXTURES = Path(__file__).parents[1] / 'real_data/fixtures'
EXPECTED = Path(__file__).parents[1] / 'real_data/expected/thiolation_plate_yield.csv'
REACTION = (
    '[c:1]1([a:2][a:3][a:4][a:5]1)[Cl,Br:6].[SH1:7][#6:8]>>'
    '[c:1]1([a:2][a:3][a:4][a:5]1)[S:7][#6:8]')


@pytest.fixture
def thiolation_layout():
    return Reaction_Array(
        FIXTURES / 'thiolation_plate_layout.csv', Transformation(REACTION),
        FIXTURES / 'thiolation_meta_data.json')


@pytest.fixture
def thiolation_results():
    expected = pd.read_csv(EXPECTED, index_col=0)
    dtype = np.dtype([
        ('quantity', float), ('rt_ms', float), ('rt_fid', float), ('flags', int)])
    results = np.zeros((8, 12), dtype=dtype)
    results['quantity'] = np.nan
    results['rt_ms'] = np.nan
    results['rt_fid'] = np.nan
    for position, row in expected.iterrows():
        index = ord(position[0]) - ord('A'), int(position[1:]) - 1
        results[index] = (
            row['Yield [%]'], row['RT-MS [min]'], row['RT-FID [min]'], 0)
    return results


def test_ord_export_builds_valid_dataset_from_study_results(
        thiolation_layout, thiolation_results, tmp_path):
    output = tmp_path / 'thiolation.pbtxt'

    dataset = Reaction_Parser.build_dataset(
        thiolation_layout, np.nan_to_num(thiolation_results['quantity']), output)

    assert dataset.name == 'pyGecko reaction array'
    assert dataset.description == 'Combinatorial reaction array exported by pyGecko.'
    assert len(dataset.reactions) == 96
    assert dataset.reactions[0].identifiers[0].value == 'A1'
    assert dataset.reactions[0].outcomes[0].products[0].measurements[0].percentage.value == 75
    assert not dataset.reactions[6].outcomes[0].products
    assert output.read_text().startswith('name: "pyGecko reaction array"')


def test_pdf_report_writes_study_results_without_package_artifacts(
        thiolation_layout, thiolation_results, tmp_path):
    output = tmp_path / 'thiolation.pdf'
    package_heatmap = Path(reports.__file__).parent / 'tmp/heatmap.png'

    PDF_Report(output, 'thiolation', thiolation_layout, thiolation_results)

    assert output.read_bytes().startswith(b'%PDF')
    assert output.stat().st_size > 1000
    assert not package_heatmap.exists()
