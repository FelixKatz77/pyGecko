'''Tests for the Injection._check_for_missing_signal diagnostic.

find_empty_ranges reports every genuine gap now that its mask is per-point, and a real MS TIC
carries dozens of them per injection. The diagnostic therefore summarises rather than itemising:
one line per injection, so a plate run stays readable.
'''

import numpy as np
import pytest

from pygecko.gc_tools.injection.fid_injection import FID_Injection


def make_injection(signal, dt=0.01, sample_name='A3'):
    '''Builds an FID_Injection whose chromatogram is then replaced by the signal under test.'''
    clean = np.array([np.arange(len(signal)) * dt, np.ones(len(signal))])
    injection = FID_Injection({'SampleName': sample_name}, clean, dt * 2)
    injection.chromatogram = np.array([np.arange(len(signal)) * dt, np.asarray(signal, float)])
    return injection


class TestCheckForMissingSignal:

    def test_silent_when_signal_is_complete(self, capsys):
        injection = make_injection([10.0] * 50)
        capsys.readouterr()
        injection._check_for_missing_signal()
        assert capsys.readouterr().out == ''

    def test_reports_one_line_per_injection(self, capsys):
        signal = [10.0] * 100
        for i in (10, 20, 21, 22, 40, 60, 61, 80):
            signal[i] = 0.0
        injection = make_injection(signal)
        capsys.readouterr()
        injection._check_for_missing_signal()
        out = capsys.readouterr().out.strip()
        assert out.count('\n') == 0, f'expected a single line, got:\n{out}'

    def test_line_carries_sample_name_count_and_worst_gap(self, capsys):
        signal = [10.0] * 100
        for i in (10, 20, 21, 22, 40):
            signal[i] = 0.0
        injection = make_injection(signal, sample_name='B7')
        capsys.readouterr()
        injection._check_for_missing_signal()
        out = capsys.readouterr().out

        assert 'B7' in out
        assert '3' in out                 # three distinct gaps
        assert '0.20' in out or '0.2' in out   # longest gap spans 0.20-0.22

    def test_counts_distinct_gaps_not_points(self, capsys):
        '''Three contiguous zero points are one gap, not three.'''
        signal = [10.0] * 50
        signal[10:13] = [0.0, 0.0, 0.0]
        injection = make_injection(signal)
        capsys.readouterr()
        injection._check_for_missing_signal()
        assert '1 ' in capsys.readouterr().out
