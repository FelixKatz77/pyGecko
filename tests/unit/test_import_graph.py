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

    def test_parsers_import_without_psims_installed(self, import_probe):
        '''psims is an optional extra, so the package must import and fail only at call time.

        The mzML writer defers its psims import, which is what lets pyGecko be installed without
        the [mzml] extra. Blocking psims in a fresh interpreter is the only way to check that,
        since pyteomics imports psims itself when it happens to be present.
        '''
        assert import_probe(
            "import sys\n"
            "class Block:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name == 'psims' or name.startswith('psims.'):\n"
            "            raise ImportError('blocked')\n"
            "        return None\n"
            "sys.meta_path.insert(0, Block())\n"
            "from pygecko.parsers.file_writers import _load_mzml_writer\n"
            "try:\n"
            "    _load_mzml_writer()\n"
            "except ImportError as error:\n"
            "    print('pyGecko[mzml]' in str(error))\n"
        ) == 'True'
