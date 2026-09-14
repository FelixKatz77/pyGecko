import subprocess

import pytest

from pygecko.parsers import msconvert_wraper


@pytest.fixture
def run_calls(monkeypatch):
    calls = []
    monkeypatch.setattr(subprocess, 'run', lambda cmd, **kwargs: calls.append(cmd))
    return calls


@pytest.fixture
def no_msconvert_on_path(monkeypatch):
    monkeypatch.delenv('PYGECKO_MSCONVERT', raising=False)
    monkeypatch.setattr(msconvert_wraper.shutil, 'which', lambda name: None)


def test_env_var_selects_executable(monkeypatch, tmp_path, run_calls, no_msconvert_on_path):
    exe = tmp_path / 'msconvert.exe'
    exe.touch()
    monkeypatch.setenv('PYGECKO_MSCONVERT', str(exe))

    msconvert_wraper.msconvert(['a.D'], str(tmp_path))

    assert run_calls == [[str(exe), '--mzML', '-o', str(tmp_path), 'a.D']]


def test_path_lookup_used_when_env_var_unset(monkeypatch, tmp_path, run_calls, no_msconvert_on_path):
    exe = tmp_path / 'msconvert'
    exe.touch()
    monkeypatch.setattr(msconvert_wraper.shutil, 'which', lambda name: str(exe) if name == 'msconvert' else None)

    msconvert_wraper.msconvert(['a.D'], str(tmp_path))

    assert run_calls[0][0] == str(exe)


def test_env_var_wins_over_path(monkeypatch, tmp_path, run_calls, no_msconvert_on_path):
    env_exe = tmp_path / 'env_msconvert'
    env_exe.touch()
    path_exe = tmp_path / 'path_msconvert'
    path_exe.touch()
    monkeypatch.setenv('PYGECKO_MSCONVERT', str(env_exe))
    monkeypatch.setattr(msconvert_wraper.shutil, 'which', lambda name: str(path_exe))

    msconvert_wraper.msconvert(['a.D'], str(tmp_path))

    assert run_calls[0][0] == str(env_exe)


def test_missing_executable_does_not_run(tmp_path, run_calls, no_msconvert_on_path, capsys):
    msconvert_wraper.msconvert(['a.D'], str(tmp_path))

    assert run_calls == []
    assert 'not found' in capsys.readouterr().out


def test_env_var_pointing_nowhere_does_not_run(monkeypatch, tmp_path, run_calls, no_msconvert_on_path):
    monkeypatch.setenv('PYGECKO_MSCONVERT', str(tmp_path / 'does_not_exist'))

    msconvert_wraper.msconvert(['a.D'], str(tmp_path))

    assert run_calls == []


def test_resolution_happens_at_call_time(monkeypatch, tmp_path, run_calls, no_msconvert_on_path):
    exe = tmp_path / 'msconvert.exe'
    exe.touch()
    msconvert_wraper.msconvert(['a.D'], str(tmp_path))
    assert run_calls == []

    monkeypatch.setenv('PYGECKO_MSCONVERT', str(exe))
    msconvert_wraper.msconvert(['a.D'], str(tmp_path))

    assert len(run_calls) == 1
