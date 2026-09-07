'''Tests for the import cost of pygecko.reaction: ORD export must stay opt-in.'''

import importlib.util

import pytest

ORD_PROBE = "import sys; {statement}; print('ord_schema' in sys.modules)"


class TestReactionPackageImports:

    def test_importing_the_package_does_not_import_ord_schema(self, import_probe):
        assert import_probe(ORD_PROBE.format(statement='import pygecko.reaction')) == 'False'

    def test_importing_the_array_module_does_not_import_ord_schema(self, import_probe):
        '''test_product_array and test_quantify_plate only need reaction.array.'''
        assert import_probe(
            ORD_PROBE.format(statement='from pygecko.reaction.array import Product_Array')
        ) == 'False'

    @pytest.mark.skipif(importlib.util.find_spec('ord_schema') is None,
                        reason='ORD export is opt-in: pip install pyGecko[ord]')
    def test_reaction_parser_is_still_importable_from_the_package(self, import_probe):
        '''The lazy attribute must still resolve where the extra is installed.'''
        assert import_probe(
            ORD_PROBE.format(statement='from pygecko.reaction import Reaction_Parser')
        ) == 'True'
