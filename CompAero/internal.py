"""Contains helper functions used by the CompAero Project"""

from math import isnan
from typing import Union

from CompAero.types import FlowState, ShockType

# Settings for printing out outputs from classes
TOTAL_WIDTH = 70  # Width of the area including  | |
INTERNAL_VALUE_WIDTH = TOTAL_WIDTH - 2  # Width of area excluding | |


def to_string(
    name: str, value: Union[str, float, int, bool], precision: int, dot_line: bool = False
) -> str:
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
    val_str = str(value)
    if isinstance(value, (int, float)):
        val_str = str(round(value, precision))

    name = name + ":"
    sep = "-" if dot_line else ""
    width = INTERNAL_VALUE_WIDTH - len(val_str)
    return f"|{name:{sep}<{width}}{val_str}|\n"


def named_subheader(name: str) -> str:
    """This generates a field which has a name in it with similiar format as to to_string()
       To be used in to seperate sub fields

    Args:
        name (str): name to print, is centered

    Returns:
        str: A formatted string with new line character on the end
    """
    name = f" {name} "
    w = INTERNAL_VALUE_WIDTH
    return f"|{name:=^{w}}|\n"


def seperator() -> str:
    """Generates and empty seperation field

    Returns:
        str: A formatted string that acts as a blank line for tables
    """
    return f"|{'':{INTERNAL_VALUE_WIDTH}}|\n"


def named_header(name: str, value: Union[float, int], precision: int) -> str:
    """Generates a title header for the table

    Args:
        name (str): name to be put in header
        value (Union[float, int, bool]): Value that is to be printed
        precision (int, optional): precision to round value to. Defaults to 4.
    Returns:
        str: A formatted string that can be used as a header
    """
    rv = round(value, precision)
    data = f" {name}: {rv:.{precision}f} "
    data_line = f"|{data:^{INTERNAL_VALUE_WIDTH}}|\n"
    f = footer()
    return f"{f}{data_line}{f}"


def footer() -> str:
    """Generates a formatted footer for the end of a table

    Returns:
        str: formatted footer
    """
    return f"|{'':=^{INTERNAL_VALUE_WIDTH}}|\n"


def check_value(*args: Union[int, float, ShockType, FlowState]) -> bool:
    """Checks to see if value is non NAN and greater than zero"""

    for arg in args:
        if isinstance(arg, (float, int)):
            if isnan(arg):
                return False

            if arg < 0.0:
                return False

        elif isinstance(arg, (FlowState, ShockType)):
            continue
        else:
            raise TypeError(
                f"Error: cannot validate variable of type {type(arg)}, \
                must be of type int, float, ShockType, or FlowState"
            )

    return True


class InvalidOptionCombinationError(Exception):
    """Thrown when the options supplied to a class are incorrect or invalid (nan)"""

    def __init__(self, *args: object) -> None:
        super().__init__(
            "Either the options suppled were the incorrect combination or not enough \
            valid arguments were supplied",
            *args,
        )


class GammaNotDefinedError(Exception):
    """Thrown when gamma is found not to be defined"""

    def __init__(self, *args: object) -> None:
        super().__init__("Gamma must be defined for the determination of flow states", *args)
