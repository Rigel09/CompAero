from enum import Enum, auto
from math import isnan
from typing import List, Union
import sys

from colorama.ansi import Style

# Settings for printing out outputs from classes
TOTAL_WIDTH = 70  # Width of the area including  | |
INTERNAL_VALUE_WIDTH = TOTAL_WIDTH - 2  # Width of area excluding | |


def value_to_string(name: str, value: Union[float, int, bool], precision: int, dot_line: bool = False) -> str:
    """This generates a professional easy to read string for a data value

    Args:
        name (str): Name that is to be printed with value
        value (Union[float, int, bool]): Value that is to be printed
        precision (int, optional): precision to round value to. Defaults to 4.
        dot_line (bool, optional): prints a dotted line between the name and the value This is used
                                    when multiple values are printed in a column. Defaults to False.

    Returns:
        str: A formatted string with new line character on the end
    """
    valString = str(round(value, precision)) if not isinstance(value, (bool, str,)) else str(value)
    name = name + ":"
    sep = "-" if dot_line else ""
    return "|{:{sep}<{width}}{}|{}\n".format(
        name, valString, Style.RESET_ALL, width=INTERNAL_VALUE_WIDTH - len(valString), sep=sep
    )


def named_subheader(name: str) -> str:
    """This generates a field which has a name in it with similiar format as to value_to_string()
       To be used in to seperate sub fields

    Args:
        name (str): name to print, is centered

    Returns:
        str: A formatted string with new line character on the end
    """
    return "|{:=^{width}}|\n".format(" " + name + " ", width=INTERNAL_VALUE_WIDTH)


def seperator() -> str:
    """Generates and empty seperation field

    Returns:
        str: A formatted string that acts as a blank line for tables
    """
    return "|{:{width}}|\n".format("", width=INTERNAL_VALUE_WIDTH)


def named_header(name: str, value: Union[float, int], precision: int) -> str:
    """Generates a title header for the table

    Args:
        name (str): name to be put in header
        value (Union[float, int, bool]): Value that is to be printed
        precision (int, optional): precision to round value to. Defaults to 4.
    Returns:
        str: A formatted string that can be used as a header
    """
    return (
        footer()
        + "|{:^{width}}|\n".format(
            " {}: {:.{precision}f} ".format(name, round(value, precision), precision=precision),
            width=INTERNAL_VALUE_WIDTH,
        )
        + footer()
    )


def footer() -> str:
    """Generates a formatted footer for the end of a table

    Returns:
        str: formatted footer
    """
    return "|{:=^{width}}|\n".format("", width=INTERNAL_VALUE_WIDTH)


def checkValue(value: Union[Union[float, int], List[Union[float, int]]]) -> bool:
    """ Checks to see if value is non NAN and greater than zero"""
    checkVal = True
    if isinstance(value, (float, int,)):
        value = [value]

    if isinstance(value, list):
        for val in value:
            checkVal = checkVal and not isnan(val)
            checkVal = checkVal and val > 0.0
    else:
        raise TypeError("{} Expected a list".format(sys._getframe().f_code.co_name))

    return checkVal


class FlowState(Enum):
    SUB_SONIC = 1
    SUPER_SONIC = 2


class ShockType(Enum):
    WEAK = auto()
    STRONG = auto()


class InvalidOptionCombinationError(Exception):
    """ Thrown when the options supplied to a class are incorrect or invalid (nan) """

    def __init__(self, *args: object) -> None:
        super().__init__(
            "Either the options suppled were the incorrect combination or not enough valid arguments were supplied",
            *args,
        )
