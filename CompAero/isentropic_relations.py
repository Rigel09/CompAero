"""This module contains the class that represents the Isentropic Flow State"""

from dataclasses import dataclass
from math import nan, sqrt

from colorama.ansi import Fore
from scipy.optimize import brenth  # type: ignore

from CompAero.greek_letters import LowerCaseGreek as lcg
from CompAero.internal import (
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    check_value,
    footer,
    named_header,
    seperator,
    to_string,
)
from CompAero.types import FlowState

ISENTROPIC_VALID_OPTIONS = [
    "gamma, mach",
    " gamma, P0/P",
    " gamma, T0/T",
    " gamma, Rho0/Rho",
    " gamma, A/A*",
]


@dataclass
class ISENTROPIC_CHOICE:
    MACH = "gamma, mach"
    P0_P = " gamma, P0/P"
    T0_T = " gamma, T0/T"
    RHO0_RHO = " gamma, Rho0/Rho"
    a_a_star = " gamma, A/A*"


class IsentropicRelations:
    # pylint: disable=too-many-instance-attributes
    """This class is a collective name space for basic calculations regarding isentropic flows.
    The constructor of this class can also determine the entire state of the flow given a partial
    state of the flow

    Args:
        gamma (float): ratio of specific heats
        mach (float, optional): Mach number of the flow. Defaults to nan.
        p0_p (float, optional): Ratio of total pressure to static pressure. Defaults to nan.
        t0_t (float, optional): Ratio of total temperature to static temperature. Defaults to nan.
        rho0_rho (float, optional): Ratio of total density to static density. Defaults to nan.
        a_a_star (float, optional): Ratio of Nozzle Area to sonic throat area. Defaults to nan.
        flow_type (FlowState, optional): States wether the flow is subsonic or supersonic.
                                        Used for Area Ratio Calculations. Defaults to
                                        FlowState.SUPER_SONIC.

    Raises:
        InvalidOptionCombinationError: Raised if invalid combination of values are given to the
                                        constructor

    Useage:
        To use this class pass gamma and one of the known parameters of the flow and the rest are
        calculated#.

    Valid_Combinations_of_Parameters:
        1: gamma, mach\n
        2: gamma, p0/p\n
        3: gamma, t0/t\n
        4: gamma, rho0/rho\n
        5: gamma, a/a*, flowtype (flow type defaults to super sonic)\n
    """

    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        p0_p: float = nan,
        t0_t: float = nan,
        rho0_rho: float = nan,
        a_a_star: float = nan,
        flow_type: FlowState = FlowState.SUPER_SONIC,
    ) -> None:
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-statements
        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ mach number of the flow """
        self.p0_p = p0_p
        """ Ratio of Total pressure to static pressure """
        self.t0_t = t0_t
        """ Ratio of total temperature to static temperature """
        self.rho0_rho = rho0_rho
        """ Ratio of total density to static density """
        self.a_a_star = a_a_star
        """ Ratio of Nozzle area to sonic nozzle area """
        self.flow_type = flow_type
        """ Type of flow which is either subsonic or supersonic (Type: flowstate) """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        # Calculate parameters based on what was passed in
        if not check_value(self.gamma):
            raise GammaNotDefinedError()

        if check_value(self.mach):
            pass
        elif check_value(self.p0_p):
            self.mach = IsentropicRelations.calc_mach_from_p0_p(self.p0_p, self.gamma)
        elif check_value(self.t0_t):
            self.mach = IsentropicRelations.calc_mach_from_t0_t(self.t0_t, self.gamma)
        elif check_value(self.rho0_rho):
            self.mach = IsentropicRelations.calc_mach_from_rho0_rho(self.rho0_rho, self.gamma)
        elif check_value(self.a_a_star):
            self.mach = IsentropicRelations.calc_mach_from_a_a_star(
                self.a_a_star, self.gamma, self.flow_type
            )
        else:
            raise InvalidOptionCombinationError()

        if check_value(self.mach):
            self.__calc_state_from_mach()

    def __calc_state_from_mach(self) -> None:
        self.t0_t = IsentropicRelations.calc_t0_t(self.mach, self.gamma)
        self.p0_p = IsentropicRelations.calc_p0_p(self.mach, self.gamma)
        self.rho0_rho = IsentropicRelations.calc_rho0_rho(self.mach, self.gamma)
        self.a_a_star = IsentropicRelations.calc_a_a_star(self.mach, self.gamma)
        self.flow_type = FlowState.SUPER_SONIC if self.mach > 1.0 else FlowState.SUB_SONIC

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Isentropic Flow State at Mach", self.mach, self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision, dot_line=True),
                to_string("p0/p", self.p0_p, self.precision),
                to_string("T0/T", self.t0_t, self.precision, dot_line=True),
                to_string(f"{lcg.rho}0/{lcg.rho}", self.rho0_rho, self.precision),
                to_string("A/A*", self.a_a_star, self.precision, dot_line=True),
                to_string("Flow Type", self.flow_type.name, self.precision),
                footer(),
            ]
        )

    @staticmethod
    def calc_t0_t(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates the ratio of Total Temperature to static temperature (T0/T)

        Args:
            mach (float): Mach number of the flow
            gamma (float): Ratio of specific heats of the flow
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: T0/T
        """

        return 1 + (gamma - 1) / 2 * pow(mach, 2) - offset

    @staticmethod
    def calc_mach_from_t0_t(t0_t: float, gamma: float) -> float:
        """Calculates the Mach number for a flow given the ratio of total tempto static temp

        Args:
            t0_t (float): Ratio of total temperature to static temperature
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return sqrt((t0_t - 1) * 2 / (gamma - 1))

    @staticmethod
    def calc_p0_p(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates the ratio of Total pressure to static pressure (P0/P)

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: P0/P
        """

        return pow((1 + (gamma - 1) / 2 * pow(mach, 2)), gamma / (gamma - 1)) - offset

    @staticmethod
    def calc_mach_from_p0_p(p0_p: float, gamma: float) -> float:
        """Calculates the Mach # for a flow given the ratio of total pressure to static pressure

        Args:
            p0_p (float): Ratio of total pressure to static pressure
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """

        return brenth(IsentropicRelations.calc_p0_p, 0, 30, args=(gamma, p0_p))  # type: ignore

    @staticmethod
    def calc_rho0_rho(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates the ratio of Total density to static density (rho0/rho)

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: rho0/rho
        """

        return pow((1 + (gamma - 1) / 2 * pow(mach, 2)), 1 / (gamma - 1)) - offset

    @staticmethod
    def calc_mach_from_rho0_rho(rho0_rho: float, gamma: float) -> float:
        """Calculates the Mach number for a flow given the ratio of total density to static density

        Args:
            rho0_rho (float): Ratio of total density to static density
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """

        return brenth(
            IsentropicRelations.calc_rho0_rho, 0, 30, args=(gamma, rho0_rho)
        )  # type: ignore

    @staticmethod
    def calc_a_a_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates the ratio of Nozzle Area to Sonic Throat Area (A/A*)

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: A/A*
        """

        gm1 = gamma - 1
        gp1 = gamma + 1
        m_sqr = pow(mach, 2)

        non_raised = 2 / gp1 * (1 + gm1 / 2 * m_sqr)

        return sqrt(pow(non_raised, gp1 / gm1) / m_sqr) - offset

    @staticmethod
    def calc_mach_from_a_a_star(
        a_a_star: float, gamma: float, flow_type: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach number for a flow given the nozzle area ratio and flow type

        Args:
            a_a_star (float): Ratio of nozzle area ratio to sonic area ratio
            gamma (float): ratio of specific heats
            flow_type (FlowState, optional): Type of flow whether it is super sonic of subsonic.
                                                Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        assert isinstance(flow_type, FlowState)
        if a_a_star == 1.0:
            return a_a_star

        if flow_type == FlowState.SUPER_SONIC:
            return brenth(
                IsentropicRelations.calc_a_a_star, 1, 30, args=(gamma, a_a_star)
            )  # type: ignore

        if flow_type == FlowState.SUB_SONIC:
            return brenth(
                IsentropicRelations.calc_a_a_star, 0.001, 1, args=(gamma, a_a_star)
            )  # type: ignore

        err = f"{Fore.RED}Unsupported flow type passed to A/A* calculations. Got {flow_type}"
        raise ValueError(err)
