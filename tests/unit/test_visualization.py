from pathlib import Path

import numpy as np

from pygecko.visualization import Visualization

from .conftest import make_fid_injection, make_peak


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


def test_chromatogram_views_write_png(tmp_path):
    first = make_fid_injection(sample_name='plate-A1', solvent_delay=0, points=500)
    second = make_fid_injection(sample_name='plate-A2', solvent_delay=0, points=500)
    for injection in (first, second):
        injection.processed_chromatogram = injection.chromatogram
        peak = make_peak(4.0)
        injection.peaks = {peak.rt: peak}

    chromatogram = tmp_path / 'chromatogram.png'
    stack = tmp_path / 'stack.png'
    Visualization.view_chromatogram(
        first, chromatogram, xlim=[3.5, 4.5], ylim=[0, None], linewidth=0.5)
    Visualization.stack_chromatograms(
        [first, second], stack, color=['#005573', '#e04214'],
        xlim=[3.5, 4.5], ylim=[0, None])

    assert chromatogram.read_bytes().startswith(b'\x89PNG')
    assert stack.read_bytes().startswith(b'\x89PNG')


def test_mass_spectrum_views_write_png(ms_peak_factory, tmp_path):
    observed = ms_peak_factory(
        5.9, {43: 25, 77: 50, 105: 100, 157: 70, 188: 20})
    reference = ms_peak_factory(
        5.9, {43: 20, 77: 55, 105: 90, 157: 75, 188: 25})
    spectrum = tmp_path / 'spectrum.png'
    comparison = tmp_path / 'comparison.png'

    Visualization.view_mass_spectrum(
        observed, spectrum, xlim=[35, 200], ylim=[0, 110], alpha=0.8)
    Visualization.compare_mass_spectra(
        (observed, reference), comparison, xlim=[35, 200], ylim=[-110, 110])

    assert spectrum.read_bytes().startswith(b'\x89PNG')
    assert comparison.read_bytes().startswith(b'\x89PNG')
