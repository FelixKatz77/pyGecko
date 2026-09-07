'''Shared fixture-path anchor for the integration suite.

The fixture paths used to be relative literals, which only resolved when pytest was invoked from
inside this directory. Anchoring them on __file__ makes the suite run identically from the repo
root, which is what [tool.pytest.ini_options] testpaths = ["tests"] implies.
'''

from pathlib import Path

FIXTURES = Path(__file__).parent / 'fixtures'


def fixture_path(*parts) -> str:
    '''Returns an absolute path inside the fixtures directory, as a string.

    The parsers take string paths, so this returns str rather than Path to keep the call sites
    identical to what they passed before.
    '''

    return str(FIXTURES.joinpath(*parts))
