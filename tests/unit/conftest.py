import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from pygecko.gc_tools.injection.fid_injection import FID_Injection
from pygecko.gc_tools.injection.injection import Injection
from pygecko.gc_tools.injection.ms_injection import MS_Injection
from pygecko.gc_tools.peak import MS_Peak
from pygecko.gc_tools.peak.peak import Peak

MS_DTYPE = dict(names=['mz', 'intensity', 'rel_intensity'], formats=['f8', 'f8', 'f8'])


def make_peak(rt, height=100.0, flags=None):
    '''Builds a bare FID-style Peak at a retention time.'''
    peak = Peak(rt, height, 0.1, np.array([rt - 0.05, rt + 0.05]), area=height)
    if flags:
        peak.flags = list(flags)
    return peak


def make_injection(peaks, sample_name='SMP-A1'):
    '''Builds a bare Injection holding the given peaks, keyed by retention time.'''
    return Injection({'SampleName': sample_name}, {p.rt: p for p in peaks})


def make_fid_chromatogram(peak_rts=(4.0, 6.0), points=4000, run_time=10.0):
    '''Builds a two-row FID chromatogram with a Gaussian peak at each given retention time.'''
    time = np.linspace(0.0, run_time, points)
    intensity = np.full(points, 5.0)
    for rt in peak_rts:
        intensity = intensity + 1000.0 * np.exp(-0.5 * ((time - rt) / 0.03) ** 2)
    return np.array([time, intensity])


def make_fid_injection(sample_name='SMP-A1', solvent_delay=1.0, **kwargs):
    '''Builds an FID_Injection over a synthetic chromatogram, ready for peak picking.'''
    return FID_Injection({'SampleName': sample_name}, make_fid_chromatogram(**kwargs), solvent_delay)


def make_mass_spectrum(mz_to_intensity):
    '''Builds a structured mass-spectrum array from an {mz: absolute intensity} mapping.'''
    intensities = np.array(list(mz_to_intensity.values()), dtype=float)
    rel = intensities / intensities.max() * 100
    return np.array(
        [(mz, i, r) for mz, i, r in zip(mz_to_intensity.keys(), intensities, rel)],
        dtype=MS_DTYPE,
    )


def make_scans(rt_ms, rows):
    '''Builds a scans DataFrame in the layout the parsers produce.

    Args:
        rt_ms (list[float]): Retention times in milliseconds, one per scan.
        rows (list[dict[int, float]]): One {nominal m/z: absolute intensity} mapping per scan.
            m/z absent from a scan are zero-filled, matching the readers.

    Returns:
        pd.DataFrame: Float millisecond index named 'retention_time', ascending integer m/z
            columns, no NaN -- the layout extract_scans_from_mzml returns.
    '''

    index = pd.Index(np.asarray(rt_ms, dtype=float), name='retention_time')
    df = pd.DataFrame(rows, index=index)
    df = df.fillna(0.0).astype(float)
    return df.reindex(sorted(df.columns), axis=1)


def make_ms_injection(peaks, sample_name='SMP-A1', scans=None):
    '''Builds an MS_Injection with a minimal two-row chromatogram and the given MS peaks.'''
    chromatogram = np.array([[0.0, 0.1, 0.2, 0.3], [1.0, 1.0, 1.0, 1.0]])
    return MS_Injection(
        {'SampleName': sample_name},
        chromatogram,
        {p.rt: p for p in peaks},
        scans=scans,
    )


@pytest.fixture
def ms_peak_factory():
    '''Returns a factory building an MS_Peak with an explicit mass spectrum.'''

    def _factory(rt, mz_to_intensity, height=100.0):
        return MS_Peak(rt, height, 0.1, np.array([rt - 0.05, rt + 0.05]),
                       make_mass_spectrum(mz_to_intensity))

    return _factory


@pytest.fixture
def import_probe():
    '''Returns a callable running Python source in a clean interpreter, returning its stdout.

    Import behaviour can only be tested out-of-process: once any test in this session has
    imported pygecko, sys.modules is warm and an in-process check passes vacuously.
    '''

    def _probe(source):
        result = subprocess.run([sys.executable, '-c', source], capture_output=True, text=True)
        assert result.returncode == 0, f'`{source}` failed:\n{result.stderr}'
        return result.stdout.strip()

    return _probe
