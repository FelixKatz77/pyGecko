from pathlib import Path

import numpy as np

from pygecko.visualization import Visualization


def test_visualize_plate_writes_structured_analysis_results(tmp_path):
    dtype = np.dtype([('quantity', float), ('flags', int)])
    results = np.zeros((2, 3), dtype=dtype)
    results['quantity'] = [[75, np.nan, 12], [0, 101, 50]]
    results['flags'][0, 2] = 1
    output = tmp_path / 'plate.png'

    Visualization.visualize_plate(
        results, output, show_flags=True,
        row_labels=['A', 'B'], col_labels=['1', '2', '3'],
        cbar_label='Conversion [%]')

    assert output.read_bytes().startswith(b'\x89PNG')
