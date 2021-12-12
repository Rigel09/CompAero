from math import isnan
from typing import List, Union


def checkValue(value: Union[float, List[float]]) -> bool:
    """ Checks to see if value is non NAN and greater than zero"""
    checkVal = True
    if isinstance(value, float):
        value = [value]

    if isinstance(value, list):
        for val in value:
            checkVal = checkVal and not isnan(val)
            checkVal = checkVal and val > 0.0

    return checkVal
