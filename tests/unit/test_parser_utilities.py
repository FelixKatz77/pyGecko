'''Tests for the parsers' filesystem helpers.'''

import pytest

from pygecko.parsers.utilities import find_directories_with_extension


def test_a_path_that_is_not_a_directory_raises(tmp_path):
    with pytest.raises(NotADirectoryError):
        find_directories_with_extension(tmp_path / 'missing', '.D')


def test_matching_directories_are_returned(tmp_path):
    (tmp_path / 'A1.D').mkdir()
    (tmp_path / 'notes.D').touch()
    assert find_directories_with_extension(tmp_path, '.D') == [tmp_path / 'A1.D']
