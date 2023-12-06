from pathlib import Path
import json

def get_num_substrates(rxn_smarts:str) -> int:

    '''
    Returns the number of substrates in a reaction SMARTS string.

    Args:
        rxn_smarts (str): Reaction SMARTS string.

    Returns:
        int: Number of substrates in the reaction.
    '''

    substrates = rxn_smarts.split('>>')[0]
    subst_list = substrates.split('.')
    num_subst = len(subst_list)
    return num_subst

def read_json(filename:Path|str) -> dict:

    '''Takes a filename of a JSON file and returns the contents of the file as a dictionary.'''

    with open(filename) as file:
        return json.load(file)
