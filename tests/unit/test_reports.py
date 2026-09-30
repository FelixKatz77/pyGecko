'''PDF_Report molecule rendering skips what Indigo cannot parse, with a warning.'''

import pytest

from pygecko.data_handling import reports
from pygecko.data_handling.reports import PDF_Report


def test_an_unparsable_smiles_warns_and_the_rest_are_rendered():
    # The report sets the output format while it builds its tables; set here, the grid can render.
    reports.indigo.setOption('render-output-format', 'png')
    report = object.__new__(PDF_Report)
    with pytest.warns(UserWarning, match='C1CC could not be rendered'):
        image = report._PDF_Report__create_common_molecules_image(['C1CC', 'CCO'], 0.5, 0.5)
    assert image is not None
