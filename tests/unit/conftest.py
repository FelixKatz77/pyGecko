import numpy as np
import pytest

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


def make_mass_spectrum(mz_to_intensity):
    '''Builds a structured mass-spectrum array from an {mz: absolute intensity} mapping.'''
    intensities = np.array(list(mz_to_intensity.values()), dtype=float)
    rel = intensities / intensities.max() * 100
    return np.array(
        [(mz, i, r) for mz, i, r in zip(mz_to_intensity.keys(), intensities, rel)],
        dtype=MS_DTYPE,
    )


def make_ms_injection(peaks, sample_name='SMP-A1'):
    '''Builds an MS_Injection with a minimal two-row chromatogram and the given MS peaks.'''
    chromatogram = np.array([[0.0, 0.1, 0.2, 0.3], [1.0, 1.0, 1.0, 1.0]])
    return MS_Injection(
        {'SampleName': sample_name},
        chromatogram,
        {p.rt: p for p in peaks},
        scans=None,
    )


@pytest.fixture
def ms_peak_factory():
    '''Returns a factory building an MS_Peak with an explicit mass spectrum.'''

    def _factory(rt, mz_to_intensity, height=100.0):
        return MS_Peak(rt, height, 0.1, np.array([rt - 0.05, rt + 0.05]),
                       make_mass_spectrum(mz_to_intensity))

    return _factory
