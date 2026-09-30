'''Tests for Transformation, the reaction-SMARTS product generator.'''

import pytest

from pygecko.reaction.transformation import Transformation


def test_a_substrate_count_mismatch_warns_and_returns_no_product():
    transformation = Transformation('[C:1][OH:2]>>[C:1][Cl]')
    with pytest.warns(UserWarning, match='Number of Substrates'):
        assert transformation(['CO', 'CCO']) is None


def test_a_matching_substrate_count_returns_the_product():
    assert Transformation('[C:1][OH:2]>>[C:1][Cl]')(['CO']) == 'CCl'
