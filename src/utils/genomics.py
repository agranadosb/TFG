from typing import Union

NUCLETODIES: list = ["A", "C", "G", "T"]

def generate_dict_values(lst: Union[tuple, list]) -> dict:
    """Function that creates a dcitionary that represents a map between a list of
    symbols and the A, C, G, T nucleotides

    For example, if our list is:
        
        ['q', 'w', 'e', 'r']
    
    the function will return:
    
        {
            'A': 'q'
            'C': 'w'
            'G': 'e'
            'T': 'r'
        }

    Parameters
    ----------
    lst: list, tuple
        List of symbols
    
    Returns
    -------
    Dictionary with the mapping
    """
    return {key: value for (key, value) in zip(NUCLETODIES, lst)}
