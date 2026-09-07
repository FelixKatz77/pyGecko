'''Tests for the import cost of pygecko.reaction: ORD export must stay opt-in.'''

import subprocess
import sys


def imports_ord_schema(statement):
    '''Runs an import in a clean interpreter, reports whether ord_schema came with it.

    A subprocess is required: once any other test has imported ord_schema, sys.modules
    is polluted for the rest of the session and the check always passes.
    '''
    code = f'import sys; {statement}; print("ord_schema" in sys.modules)'
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                            check=True)
    return result.stdout.strip().splitlines()[-1] == 'True'


class TestReactionPackageImports:

    def test_importing_the_package_does_not_import_ord_schema(self):
        assert not imports_ord_schema('import pygecko.reaction')

    def test_importing_the_array_module_does_not_import_ord_schema(self):
        '''test_product_array and test_quantify_plate only need reaction.array.'''
        assert not imports_ord_schema('from pygecko.reaction.array import Product_Array')

    def test_reaction_parser_is_still_importable_from_the_package(self):
        assert imports_ord_schema('from pygecko.reaction import Reaction_Parser')
