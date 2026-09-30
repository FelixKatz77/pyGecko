'''MS_Base_Parser file discovery: which entries of a raw directory become injections.'''

from types import SimpleNamespace

from pygecko.parsers.ms_base_parser import MS_Base_Parser


def record_loads(monkeypatch):
    loaded = []

    def fake_load_injection(raw_data_path, pos=False, temp_dir=None):
        loaded.append(raw_data_path.name)
        return SimpleNamespace(sample_name=raw_data_path.stem)

    monkeypatch.setattr(MS_Base_Parser, 'load_injection', staticmethod(fake_load_injection))
    return loaded


def test_a_name_merely_ending_in_cdf_is_not_loaded(tmp_path, monkeypatch):
    # The pattern '*cdf' (no dot) matched x_notcdf, which then went to msconvert.
    (tmp_path / 'x_spectra.cdf').touch()
    (tmp_path / 'x_notcdf').touch()
    loaded = record_loads(monkeypatch)
    MS_Base_Parser.load_sequence(tmp_path)
    assert loaded == ['x_spectra.cdf']


def test_every_supported_format_is_found_whatever_its_case(tmp_path, monkeypatch):
    (tmp_path / 'A1.D').mkdir()
    for name in ('A2.mzml', 'A3.mzXML', 'A4_spectra.CDF', 'notes.txt'):
        (tmp_path / name).touch()
    loaded = record_loads(monkeypatch)
    MS_Base_Parser.load_sequence(tmp_path)
    assert sorted(loaded) == ['A1.D', 'A2.mzml', 'A3.mzXML', 'A4_spectra.CDF']
