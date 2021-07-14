from math import isnan
from typing import Union

def checkValue(var: float) -> bool:
    ''' Checks to see if value is non NAN and greater than zero'''
    if isnan(var): return False

    if var < 0: return False

    return True