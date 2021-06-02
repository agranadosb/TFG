from typing import Union

NUCLETODIES: list = ["A", "C", "G", "T"]
""" List of nucleotides """


def generate_dict_values(
    lst: Union[tuple, list], reference: Union[tuple, list] = NUCLETODIES
) -> dict:
    """Function that creates a dcitionary that represents a map between a list of
    symbols and the reference nucleotides (A, C, G and T for example)

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
    reference: list, tuple = ["A", "C", "G", "T"]
        List of reference symbols

    Returns
    -------
    Dictionary with the mapping
    """
    return {key: value for (key, value) in zip(reference, lst)}
