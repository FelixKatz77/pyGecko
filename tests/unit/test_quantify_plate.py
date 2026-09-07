'''Tests for the plate shape of the single-detector Analysis.quantify_plate path.'''

import numpy as np
import pytest

from pygecko.analysis import Analysis
from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.sequence.fid_sequence import FID_Sequence
from pygecko.reaction.array import Product_Array


def make_fid_injection(pos):
    '''Builds a peakless FID_Injection at a plate position.

    Peakless is deliberate: flag_peak then finds nothing, every well reports NaN, and the test
    exercises the grid assembly without dragging in RDKit-based Polyarc quantification.
    '''
    time_axis = np.arange(0, 1.0, 0.01)
    injection = FID_Injection({'SampleName': pos},
                              np.array([time_axis, np.ones_like(time_axis)]),
                              0.05, pos=True)
    injection.peaks = {}
    return injection


def make_fid_sequence(wells):
    return FID_Sequence({}, {pos: make_fid_injection(pos) for pos in wells})


@pytest.fixture
def layout_11x3(tmp_path):
    '''An 11x3 (A1..K3) layout, the shape of the split-GC reference sequence.'''
    path = tmp_path / 'products.csv'
    rows = '\n'.join(','.join(f'C{i}{c}' for c in '123') for i in range(11))
    path.write_text(f'1,2,3\n{rows}\n')
    return Product_Array(path)


class TestQuantifyPlateShape:

    def test_derives_shape_from_layout(self, layout_11x3):
        wells = [f'{row}{col}' for row in 'ABCDEFGHIJK' for col in '123']
        results = Analysis.quantify_plate(make_fid_sequence(wells), 5.0, layout=layout_11x3)
        assert results.shape == (11, 3)

    def test_keeps_every_well_of_a_non_8x12_plate(self, layout_11x3):
        wells = [f'{row}{col}' for row in 'ABCDEFGHIJK' for col in '123']
        results = Analysis.quantify_plate(make_fid_sequence(wells), 5.0, layout=layout_11x3)
        assert results['quantity'].size == len(wells)

    def test_defaults_to_8x12_without_a_layout(self):
        wells = [f'{row}{col}' for row in 'ABCDEFGH' for col in range(1, 13)]
        results = Analysis.quantify_plate(make_fid_sequence(wells), 5.0)
        assert results.shape == (8, 12)

    def test_result_dtype_is_unchanged(self, layout_11x3):
        wells = [f'{row}{col}' for row in 'ABCDEFGHIJK' for col in '123']
        results = Analysis.quantify_plate(make_fid_sequence(wells), 5.0, layout=layout_11x3)
        assert results.dtype.names == ('quantity', 'rt_fid', 'flags')
