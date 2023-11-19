"""Contains dataclass for converting Greek letters to strings"""

from dataclasses import dataclass


@dataclass
class UpperCaseGreek:
    """Class that maps uppercase unicode greek letter values"""

    # pylint: disable=too-many-instance-attributes

    alpha: str = "\u0391"
    beta: str = "\u0392"
    gamma: str = "\u0393"
    delta: str = "\u0394"
    epsilon: str = "\u0395"
    zeta: str = "\u0396"
    eta: str = "\u0397"
    theta: str = "\u0398"
    iota: str = "\u0399"
    kappa: str = "\u039A"
    lamda: str = "\u039B"
    mu: str = "\u039C"
    nu: str = "\u039D"
    xi: str = "\u039E"
    omicron: str = "\u039F"
    pi: str = "\u03A0"
    rho: str = "\u03A1"
    sigma: str = "\u03A3"
    tau: str = "\u03A4"
    upsilon: str = "\u03A5"
    phi: str = "\u03A6"
    chi: str = "\u03A7"
    psi: str = "\u03A8"
    omega: str = "\u03A9"
    symbol: str = "\u03F4"


@dataclass
class LowerCaseGreek:
    """Class that maps lowercase unicode greek letter values"""

    # pylint: disable=too-many-instance-attributes

    alpha: str = "\u03B1"
    beta: str = "\u03B2"
    gamma: str = "\u03B3"
    delta: str = "\u03B4"
    epsilon: str = "\u03B5"
    zeta: str = "\u03B6"
    eta: str = "\u03B7"
    theta: str = "\u03B8"
    iota: str = "\u03B9"
    kappa: str = "\u03BA"
    lamda: str = "\u03BB"
    mu: str = "\u03BC"
    nu: str = "\u03BD"
    xi: str = "\u03BE"
    omicron: str = "\u03BF"
    pi: str = "\u03C0"
    rho: str = "\u03C1"
    sigma: str = "\u03C2"
    tau: str = "\u03C3"
    upsilon: str = "\u03C4"
    phi: str = "\u03C5"
    chi: str = "\u03C6"
    psi: str = "\u03C7"
    omega: str = "\u03C8"
    symbol: str = "\u03C9"


@dataclass
class Misc:
    """Class that stores misc characters such as the degree symbol"""

    degree: str = "\u00b0"
