'''Analysis_Settings.pop: only None means a setting is unset.

Every setting update() has not touched is None, so None is the one value that means "use the
computed default". An explicit 0 or () is a choice the caller made and must reach the algorithm.
'''

import pytest

from pygecko.gc_tools.analysis.analysis_settings import Analysis_Settings


@pytest.mark.parametrize('key, value', [
    ('height', 0),
    ('width', 0),
    ('prominence_fid', 0),
    ('time_range', ()),
])
class TestAnExplicitFalsyValueIsKept:

    def test_pop_returns_the_configured_value(self, key, value):
        settings = Analysis_Settings()
        settings.update(**{key: value})
        assert settings.pop(key, 99) == value

    def test_resolved_records_the_configured_value(self, key, value):
        settings = Analysis_Settings()
        settings.update(**{key: value})
        settings.pop(key, 99)
        assert settings._resolved[key] == value


def test_an_unset_setting_resolves_to_the_default():
    settings = Analysis_Settings()
    assert settings.pop('height', 99) == 99
