import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.figure import Figure

from pygecko.visualization import Visualization

from .conftest import make_fid_injection, make_peak


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


@pytest.fixture
def plate_results():
    dtype = np.dtype([('quantity', float), ('flags', int)])
    results = np.zeros((2, 3), dtype=dtype)
    results['quantity'] = [[75, np.nan, 12], [0, 101, 50]]
    results['flags'][0, 2] = 1
    return results


@pytest.fixture
def injections():
    '''Two FID injections with a processed chromatogram and one peak each, ready to plot.'''
    first = make_fid_injection(sample_name='plate-A1', solvent_delay=0, points=500)
    second = make_fid_injection(sample_name='plate-A2', solvent_delay=0, points=500)
    for injection in (first, second):
        injection.chromatogram.processed = injection.chromatogram.intensity
        peak = make_peak(4.0)
        injection.peaks = {peak.rt: peak}
    return first, second


@pytest.fixture
def spectra(ms_peak_factory):
    observed = ms_peak_factory(5.9, {43: 25, 77: 50, 105: 100, 157: 70, 188: 20})
    reference = ms_peak_factory(5.9, {43: 20, 77: 55, 105: 90, 157: 75, 188: 25})
    return observed, reference


def draw_each(plate_results, injections, spectra):
    '''Returns the figure of every Visualization method, keyed by method name.'''
    return {
        'visualize_plate': Visualization.visualize_plate(
            plate_results, show_flags=True, row_labels=['A', 'B'], col_labels=['1', '2', '3'],
            cbar_label='Conversion [%]'),
        'view_chromatogram': Visualization.view_chromatogram(
            injections[0], xlim=[3.5, 4.5], ylim=[0, None], linewidth=0.5),
        'stack_chromatograms': Visualization.stack_chromatograms(
            list(injections), color=['#005573', '#e04214'], xlim=[3.5, 4.5], ylim=[0, None]),
        'view_mass_spectrum': Visualization.view_mass_spectrum(
            spectra[0], xlim=[35, 200], ylim=[0, 110], alpha=0.8),
        'compare_mass_spectra': Visualization.compare_mass_spectra(
            spectra, xlim=[35, 200], ylim=[-110, 110]),
    }


def test_every_method_returns_a_figure(plate_results, injections, spectra):
    for name, figure in draw_each(plate_results, injections, spectra).items():
        assert isinstance(figure, Figure), name


def test_every_figure_can_be_saved_by_the_caller(plate_results, injections, spectra, tmp_path):
    for name, figure in draw_each(plate_results, injections, spectra).items():
        output = tmp_path / f'{name}.png'
        figure.savefig(output)
        assert output.read_bytes().startswith(b'\x89PNG'), name


def test_closing_the_returned_figures_leaves_none_open(plate_results, injections, spectra):
    plt.close('all')
    for figure in draw_each(plate_results, injections, spectra).values():
        plt.close(figure)
    assert plt.get_fignums() == []


def test_the_style_is_applied_to_the_figure_text(injections):
    # font.size 12 instead of matplotlib's 10. Tick labels are created lazily, so they show whether
    # the style reached text built after the plotting calls, not just the calls themselves.
    ax = Visualization.view_chromatogram(injections[0]).axes[0]
    assert ax.title.get_fontsize() == pytest.approx(14.4)   # 'large' = 1.2 x font.size
    assert {label.get_fontsize() for label in ax.get_xticklabels()} == {12}


def test_importing_the_module_leaves_global_rcparams_unchanged(import_probe):
    # Out of process: once any test has imported pygecko.visualization the import is a no-op.
    out = import_probe(
        'import matplotlib.pyplot as plt; '
        'before = {k: plt.rcParams[k] for k in ("font.family", "font.sans-serif", "font.size")}; '
        'import pygecko.visualization; '
        'print(before == {k: plt.rcParams[k] for k in before})')
    assert out == 'True'


class TestStackChromatograms:

    def test_limits_apply_to_every_trace(self, injections):
        figure = Visualization.stack_chromatograms(list(injections), xlim=[3.5, 4.5], ylim=[0, 2000])
        assert [ax.get_xlim() for ax in figure.axes] == [(3.5, 4.5)] * 2
        assert [ax.get_ylim() for ax in figure.axes] == [(0, 2000)] * 2

    def test_raw_applies_to_every_trace(self, injections):
        for injection in injections:
            injection.chromatogram.processed = np.zeros_like(injection.chromatogram.intensity)
        figure = Visualization.stack_chromatograms(list(injections), raw=True)
        for ax, injection in zip(figure.axes, injections):
            assert ax.lines[0].get_ydata().max() == pytest.approx(injection.chromatogram.intensity.max())


def test_the_injection_and_peak_shortcuts_return_the_figure(injections, spectra):
    assert isinstance(injections[0].view_chromatogram(xlim=[3.5, 4.5]), Figure)
    assert isinstance(spectra[0].view_mass_spectrum(), Figure)
