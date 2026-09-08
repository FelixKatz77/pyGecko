from pygecko.reaction.layout import Combinatorial_Layout
from pygecko.reaction.transformation import Transformation
from pygecko.reaction.array import Reaction_Array, Product_Array
from pygecko.reaction.utilities import get_num_substrates, read_json

__all__ = ['Combinatorial_Layout', 'Transformation', 'Reaction_Array', 'Product_Array',
           'get_num_substrates', 'read_json', 'Reaction_Parser']


def __getattr__(name):
    '''Defers Reaction_Parser so ord_schema stays optional (pip install pyGecko[ord]).'''
    if name == 'Reaction_Parser':
        from pygecko.reaction.reaction_parser import Reaction_Parser
        return Reaction_Parser
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
