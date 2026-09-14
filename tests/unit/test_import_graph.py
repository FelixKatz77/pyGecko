'''Tests that each public subpackage imports standalone, with no import-order dependency.'''

import pytest

ENTRY_POINTS = [
    'import pygecko.visualization',
    'from pygecko.visualization import Visualization',
    'import pygecko.visualization.utilities',
    'import pygecko.data_handling',
    'from pygecko.data_handling import PDF_Report',
    'import pygecko.analysis',
    'import pygecko.gc_tools',
    'import pygecko.parsers',
    'import pygecko.reaction',
    'import pygecko.reaction; import pygecko.visualization',
    'import pygecko.data_handling; import pygecko.gc_tools',
]


class TestEntryPointsImportStandalone:

    @pytest.mark.parametrize('statement', ENTRY_POINTS)
    def test_entry_point_imports_in_a_fresh_interpreter(self, statement, import_probe):
        import_probe(statement)

    def test_importing_gc_tools_does_not_import_visualization(self, import_probe):
        '''The upward edge to Visualization must stay deferred to call time.'''
        assert import_probe(
            "import sys, pygecko.gc_tools; print('pygecko.visualization' in sys.modules)"
        ) == 'False'
