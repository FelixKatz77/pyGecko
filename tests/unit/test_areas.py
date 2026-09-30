'''FID areas are physical: intensity x minutes, independent of the sampling rate.

simpson was called without x, so an area came out in intensity x scans and doubled when the same
peak was acquired at twice the scan rate. Ratios, and so yields, were unaffected, but the reported
areas could not be compared across methods with different scan rates.
'''

import numpy as np
import pytest

from pygecko.gc_tools import Chromatogram
from pygecko.gc_tools.injection.fid_injection import FID_Injection

HEIGHT = 1000.0
SIGMA = 0.03   # minutes
RT = 6.0


def gaussian_injection(points):
    time = np.linspace(0.0, 10.0, points)
    intensity = 5.0 + HEIGHT * np.exp(-0.5 * ((time - RT) / SIGMA) ** 2)
    injection = FID_Injection({'SampleName': 'A1'}, Chromatogram(time, intensity, kind='FID'), 0.001)
    injection.pick_peaks()
    return injection


def the_area(injection):
    (peak,) = injection.peaks.values()
    return peak.area


def test_a_gaussian_integrates_to_height_times_sigma_times_root_two_pi():
    assert the_area(gaussian_injection(4000)) == pytest.approx(HEIGHT * SIGMA * np.sqrt(2 * np.pi), rel=0.01)


def test_twice_the_scan_rate_gives_the_same_area():
    assert the_area(gaussian_injection(8000)) == pytest.approx(the_area(gaussian_injection(4000)), rel=0.01)


def test_integrate_agrees_with_pick_peaks_in_the_same_units():
    injection = gaussian_injection(4000)
    picked = the_area(injection)
    injection.integrate()
    assert the_area(injection) == pytest.approx(picked, rel=1e-12)
