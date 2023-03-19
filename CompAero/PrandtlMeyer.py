from math import atan, sqrt, nan, pow, radians, degrees
from types import DynamicClassAttribute
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from enum import Enum
from CompAero.internal import (
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    checkValue,
    footer,
    named_header,
    named_subheader,
    seperator,
    to_string,
)
from CompAero.ObliqueShockRelations import ObliqueShockRelations
from CompAero.greek_letters import LowerCaseGreek as lcg, Misc

class PrandtlMeyerChoice(Enum):
    GAMMA_MACH = "gamma, mach"
    GAMMA_NU = "gamma, nu"
    GAMMA_MU = "gamma, mu"
    GAMMA_DEFLECTION_DWN_STRM_MACH = "gamma, deflection angle, dwnStrm_mach"
    GAMMA_DEFLECTION_DWN_STRM_MU = "gamma, deflection angle, dwnStrm_mu"
    GAMMA_DEFLECTION_DWN_STRM_NU = "gamma, deflection angle, dwnStrm_nu"
    
PRANDTL_MEYER_OPTIONS = [x.value for x in PrandtlMeyerChoice]

class PrandtlMeyer:
    """ This class is a collective name space for basic calculations regarding Prandtl Meyer flows. 
    The constructor of this class can also determine the entire state of the flow given a partial state of the flow 

    Args:
        gamma (float): Ratio of specific heats
        mach (float, optional): Mach number of the flow. Defaults to nan.
        nu (float, optional): Prandtl Meyer function. Defaults to nan.
        mu (float, optional): Mach Wave angle. Defaults to nan.
        dwnstreamNu (float, optional): down stream prandtl meyer function. Defaults to nan.
        dwnStreamMu (float, optional): down stream mach wave angle. Defaults to nan.
        dwnStreamMach (float, optional): down stream mach number. Defaults to nan.
        deflectionAngle (float, optional): deflection angle seen by the flow. Defaults to nan.
        inDegrees (bool, optional): True if angles passed in are in degrees. Doesnt convert output to defrees. Defaults to True.

    Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and flow state cannot be determined
        
    Useage:
    To use this class pass gamma and one of the known parameters of the flow and the rest are calculated. 

    Valid_Combinations_of_Parameters:
        1: gamma, mach\n
        2: gamma, nu\n
        3: gamma, mu\n
        4: gamma, deflection angle, dwnStrm_mach\n
        5: gamma, deflection angle, dwnStrm_mu\n
        6: gamma, deflection angle, dwnStrm_nu\n
    """

    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        nu: float = nan,
        mu: float = nan,
        dwnstreamNu: float = nan,
        dwnStreamMu: float = nan,
        dwnStreamMach: float = nan,
        deflectionAngle: float = nan,
        inDegrees: bool = True,
    ) -> None:
        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ Mach number of the flow """
        self.nu = nu
        """ The Prandtl Meyer function value """
        self.mu = mu
        """ The angle of the mach wave"""
        self.deflectionAngle = deflectionAngle
        """ Angle of the flow deflection angle"""
        self._degrees = inDegrees
        """ Angles passed to the constructor are in degrees"""
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        self.dwmStrm_nu = dwnstreamNu
        """ The Prandtl Meyer function value after the expansion wave """
        self.dwmStrm_mu = dwnStreamMu
        """ The angle of the mach wave after the expansion wave"""
        self.dwmStrm_mach = dwnStreamMach
        """ Mach number of the flow after the expansion wave """

        if not self._degrees:
            self.deflectionAngle = degrees(self.deflectionAngle)

        if not checkValue(self.gamma):
            raise GammaNotDefinedError()

        if checkValue(self.mach):
            pass
        elif checkValue(self.nu):
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        elif checkValue(self.mu):
            self.mach = ObliqueShockRelations.calc_mach_from_mach_wave_angle(radians(self.mu))

        elif checkValue(self.deflectionAngle) and checkValue(self.dwmStrm_mach):
            self.dwmStrm_nu = PrandtlMeyer.calc_nu(self.dwmStrm_mach, self.gamma)
            self.nu = self.dwmStrm_nu - self.deflectionAngle
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        elif checkValue(self.deflectionAngle) and checkValue(self.dwmStrm_nu):
            self.nu = self.dwmStrm_nu - self.deflectionAngle
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        elif checkValue(self.deflectionAngle) and checkValue(self.dwmStrm_mu):
            self.dwmStrm_mach = ObliqueShockRelations.calc_mach_from_mach_wave_angle(radians(self.dwmStrm_mu))
            self.dwmStrm_nu = PrandtlMeyer.calc_nu(self.dwmStrm_mach, self.gamma)
            self.nu = self.dwmStrm_nu - self.deflectionAngle
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        else:
            raise InvalidOptionCombinationError()

        if checkValue(self.mach) and self.mach >= 1.0:
            self.__calculateState()

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Prandtl Relations at Mach", self.mach, self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string(lcg.nu, self.nu, self.precision, dot_line=True),
                to_string(lcg.mu, self.mu, self.precision),
                seperator(),
                named_subheader("Downstream Conditions"),
                to_string("Mach", self.dwmStrm_mach, self.precision),
                to_string(lcg.nu, self.dwmStrm_nu, self.precision, dot_line=True),
                to_string(lcg.mu, self.dwmStrm_mu, self.precision),
                to_string(
                    "Flow Deflection Angle [{}]".format(lcg.theta),
                    self.deflectionAngle,
                    self.precision,
                    dot_line=True,
                ),
                footer(),
            ]
        )

    def __calculateState(self) -> None:
        self.nu = PrandtlMeyer.calc_nu(self.mach, self.gamma)
        self.mu = degrees(ObliqueShockRelations.calc_mach_wave_angle(self.mach))

        if not checkValue(self.deflectionAngle):
            return

        self.dwmStrm_nu = self.deflectionAngle + self.nu
        self.dwmStrm_mach = PrandtlMeyer.calc_mach_from_nu(self.dwmStrm_nu, self.gamma)
        self.dwmStrm_mu = degrees(ObliqueShockRelations.calc_mach_wave_angle(self.dwmStrm_mach))

    @staticmethod
    def calc_nu(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates the prandtl meyer function value (nu)

        Args:
            mach (float): mach number of the the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: nu
        """
        if mach <= 1.0:
            return 0.0

        gp1 = gamma + 1
        gm1 = gamma - 1
        mSqrMinus1 = pow(mach, 2) - 1
        return degrees(sqrt(gp1 / gm1) * atan(sqrt(gm1 / gp1 * mSqrMinus1)) - atan(sqrt(mSqrMinus1))) - offset

    @staticmethod
    def calc_mach_from_nu(nu: float, gamma: float) -> float:
        """ Calculates the mach number based on a prandtl meyer function value

        Args:
            nu (float): prandtl meyer function value
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        if nu <= 0.0:
            return 1.0

        return brenth(PrandtlMeyer.calc_nu, 1 + 1e-9, 30, args=(gamma, nu))

