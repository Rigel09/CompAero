"""
This class defines the relations for the Prandtl Meyer flows
"""

from enum import Enum
from math import atan, degrees, nan, radians, sqrt

from scipy.optimize import brenth  # type: ignore

from CompAero.greek_letters import LowerCaseGreek as lcg
from CompAero.internal import (
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    check_value,
    footer,
    named_header,
    named_subheader,
    seperator,
    to_string,
)
from CompAero.oblique_shock_relations import ObliqueShockRelations


class PrandtlMeyerChoice(Enum):
    """List valid initialization options for prandtl meyer flows"""

    GAMMA_MACH = "gamma, mach"
    GAMMA_NU = "gamma, nu"
    GAMMA_MU = "gamma, mu"
    GAMMA_DEFLECTION_DWN_STRM_MACH = "gamma, deflection angle, dwn_strm_mach"
    GAMMA_DEFLECTION_DWN_STRM_MU = "gamma, deflection angle, dwn_strm_mu"
    GAMMA_DEFLECTION_DWN_STRM_NU = "gamma, deflection angle, dwn_strm_nu"


PRANDTL_MEYER_OPTIONS = [x.value for x in PrandtlMeyerChoice]


class PrandtlMeyer:
    # pylint: disable=too-many-instance-attributes
    """
    This class is a collective name space for basic calculations regarding Prandtl Meyer flows.
    The constructor of this class determines the entire state of the flow from a partial state

    Args:
        gamma (float): Ratio of specific heats
        mach (float, optional): Mach number of the flow. Defaults to nan.
        nu (float, optional): Prandtl Meyer function. Defaults to nan.
        mu (float, optional): Mach Wave angle. Defaults to nan.
        down_stream_nu (float, optional): down stream prandtl meyer function. Defaults to nan.
        down_stream_mu (float, optional): down stream mach wave angle. Defaults to nan.
        down_stream_mach (float, optional): down stream mach number. Defaults to nan.
        deflection_angle (float, optional): deflection angle seen by the flow. Defaults to nan.
        in_degrees (bool, optional): Flags that angles passed in are in degrees. Doesnt convert
                                    output to degrees. Defaults to True.

    Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and
                                        flow state cannot be determined

    Useage:
    To use this class pass gamma and one other known parameter of the flow

    Valid_Combinations_of_Parameters:
        1: gamma, mach\n
        2: gamma, nu\n
        3: gamma, mu\n
        4: gamma, deflection angle, dwn_strm_mach\n
        5: gamma, deflection angle, dwn_strm_mu\n
        6: gamma, deflection angle, dwn_strm_nu\n
    """

    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        nu: float = nan,
        mu: float = nan,
        down_stream_nu: float = nan,
        down_stream_mu: float = nan,
        down_stream_mach: float = nan,
        deflection_angle: float = nan,
        in_degrees: bool = True,
    ) -> None:
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-statements
        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ Mach number of the flow """
        self.nu = nu
        """ The Prandtl Meyer function value """
        self.mu = mu
        """ The angle of the mach wave"""
        self.deflection_angle = deflection_angle
        """ Angle of the flow deflection angle"""
        self._degrees = in_degrees
        """ Angles passed to the constructor are in degrees"""
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        self.down_stream_nu = down_stream_nu
        """ The Prandtl Meyer function value after the expansion wave """
        self.down_stream_mu = down_stream_mu
        """ The angle of the mach wave after the expansion wave"""
        self.down_stream_mach = down_stream_mach
        """ Mach number of the flow after the expansion wave """

        if not self._degrees:
            self.deflection_angle = degrees(self.deflection_angle)

        if not check_value(self.gamma):
            raise GammaNotDefinedError()

        if check_value(self.mach):
            pass
        elif check_value(self.nu):
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        elif check_value(self.mu):
            self.mach = ObliqueShockRelations.calc_mach_from_mach_wave_angle(radians(self.mu))

        elif check_value(self.deflection_angle) and check_value(self.down_stream_mach):
            self.down_stream_nu = PrandtlMeyer.calc_nu(self.down_stream_mach, self.gamma)
            self.nu = self.down_stream_nu - self.deflection_angle
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        elif check_value(self.deflection_angle) and check_value(self.down_stream_nu):
            self.nu = self.down_stream_nu - self.deflection_angle
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        elif check_value(self.deflection_angle) and check_value(self.down_stream_mu):
            self.down_stream_mach = ObliqueShockRelations.calc_mach_from_mach_wave_angle(
                radians(self.down_stream_mu)
            )
            self.down_stream_nu = PrandtlMeyer.calc_nu(self.down_stream_mach, self.gamma)
            self.nu = self.down_stream_nu - self.deflection_angle
            self.mach = PrandtlMeyer.calc_mach_from_nu(self.nu, self.gamma)

        else:
            raise InvalidOptionCombinationError()

        if check_value(self.mach) and self.mach >= 1.0:
            self._calc_state()

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
                to_string("Mach", self.down_stream_mach, self.precision),
                to_string(lcg.nu, self.down_stream_nu, self.precision, dot_line=True),
                to_string(lcg.mu, self.down_stream_mu, self.precision),
                to_string(
                    f"Flow Deflection Angle [{lcg.theta}]",
                    self.deflection_angle,
                    self.precision,
                    dot_line=True,
                ),
                footer(),
            ]
        )

    def _calc_state(self) -> None:
        self.nu = PrandtlMeyer.calc_nu(self.mach, self.gamma)
        self.mu = degrees(ObliqueShockRelations.calc_mach_wave_angle(self.mach))

        if not check_value(self.deflection_angle):
            return

        self.down_stream_nu = self.deflection_angle + self.nu
        self.down_stream_mach = PrandtlMeyer.calc_mach_from_nu(self.down_stream_nu, self.gamma)
        self.down_stream_mu = degrees(
            ObliqueShockRelations.calc_mach_wave_angle(self.down_stream_mach)
        )

    @staticmethod
    def calc_nu(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates the prandtl meyer function value (nu)

        Args:
            mach (float): mach number of the the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding. Defaults to 0.0.

        Returns:
            float: nu
        """
        if mach <= 1.0:
            return 0.0

        gp1 = gamma + 1
        gm1 = gamma - 1
        m_sqr_minus_1 = pow(mach, 2) - 1
        return (
            degrees(
                sqrt(gp1 / gm1) * atan(sqrt(gm1 / gp1 * m_sqr_minus_1)) - atan(sqrt(m_sqr_minus_1))
            )
            - offset
        )

    @staticmethod
    def calc_mach_from_nu(nu: float, gamma: float) -> float:
        """Calculates the mach number based on a prandtl meyer function value

        Args:
            nu (float): prandtl meyer function value
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        if nu <= 0.0:
            return 1.0

        return brenth(PrandtlMeyer.calc_nu, 1 + 1e-9, 30, args=(gamma, nu))  # type: ignore
