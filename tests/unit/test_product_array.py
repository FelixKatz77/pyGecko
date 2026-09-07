'''Tests for Product_Array's product lookup API.'''

import pytest

from pygecko.reaction.array import Product_Array


@pytest.fixture
def layout_file(tmp_path):
    path = tmp_path / 'products.csv'
    path.write_text('1,2,3\nCCO,CCC,CCN\nc1ccccc1,CC=O,CCBr\n')
    return path


class TestProductArray:

    def test_get_product_returns_layout_entry(self, layout_file):
        array = Product_Array(layout_file)
        assert array.get_product('A1') == 'CCO'
        assert array.get_product('B3') == 'CCBr'

    def test_get_product_matches_getitem(self, layout_file):
        array = Product_Array(layout_file)
        for pos in ('A1', 'A2', 'A3', 'B1', 'B2', 'B3'):
            assert array.get_product(pos) == array[pos]
