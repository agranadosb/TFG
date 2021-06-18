from math import gcd
from typing import Union


def _lcm(number_1: int, number_2: int) -> int:
    """Computes the lcm between two numbers.

    Parameters
    ----------
    number_1: int
        Number 1
    number_2: int
        Number 2
    
    Returns
    -------
    Lcm
    """
    return number_1 * number_2 // gcd(number_1, number_2)

def lcm(numbers: Union[list, tuple]) -> int:
    """Computes the lcm of a list of numbers.

    Parameters
    ----------
    numbers: list, tuple
        List of numbers
    
    Returns
    -------
    Lcm
    """
    lcm_value = 1
    for number in numbers:
        lcm_value = _lcm(lcm_value, number)
    return lcm_value
