'''FID_Base_Parser file handling: which files are read and how.'''

from types import SimpleNamespace

import numpy as np
import pytest

from pygecko.parsers.fid_base_parser import FID_Base_Parser


class TestReadXyArray:

    def test_an_unsupported_suffix_raises(self, tmp_path):
        path = tmp_path / 'a.txt'
        path.write_text('0.0\t1.0\n0.1\t2.0\n')
        with pytest.raises(ValueError, match='a.txt'):
            FID_Base_Parser.read_xy_array(path)

    @pytest.mark.parametrize('name, content', [
        ('a.csv', '0.0,1.0\n0.1,2.0\n0.2,3.0\n'),
        ('a.XY', '0.0\t1.0\n0.1\t2.0\n0.2\t3.0\n'),
        ('a.CSV', '0.0,1.0\n0.1,2.0\n0.2,3.0\n'),
    ], ids=['csv', 'XY', 'CSV'])
    def test_the_suffix_is_matched_case_insensitively(self, tmp_path, name, content):
        path = tmp_path / name
        path.write_text(content)
        np.testing.assert_allclose(FID_Base_Parser.read_xy_array(path), [[0.0, 0.1, 0.2], [1.0, 2.0, 3.0]])


class TestLoadSequenceDiscovery:

    def test_every_supported_file_is_loaded_once_whatever_its_case(self, tmp_path, monkeypatch):
        for name in ('a.csv', 'b.XY', 'c.txt'):
            (tmp_path / name).touch()
        loaded = []

        def fake_load_injection(xy_file, solvent_delay, pos=False):
            loaded.append(xy_file.name)
            return SimpleNamespace(sample_name=xy_file.stem)

        monkeypatch.setattr(FID_Base_Parser, 'load_injection', staticmethod(fake_load_injection))
        FID_Base_Parser.load_sequence(tmp_path, 1.0)
        assert sorted(loaded) == ['a.csv', 'b.XY']
