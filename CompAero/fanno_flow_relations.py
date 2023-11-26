"""Provides a class around the Fanno Flow State"""

from enum import Enum
from math import log, nan, sqrt

from colorama import Fore
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
from CompAero.types import FlowState


class FannoFlowChoice(Enum):
    GAMMA_MACH = "gamma, mach"
    GAMMA_T_T_ST = "gamma, T/T*"
    GAMMA_P_P_ST = "gamma, P/P*"
    GAMMA_RHO_RHO_ST = "gamma, rho/rho*"
    GAMMA_PO_PO_ST = "gamma, P0/P0*"
    GAMMA_4FLSTD_FLOW_TYPE = "gamma, 4FL*/D, flowtype"
    GAMMA_U_U_ST = "gamma, U/U*"


FANNO_FLOW_VALID_OPTIONS = [x.value for x in FannoFlowChoice]


class FannoFlowRelations:
    # pylint: disable=too-many-instance-attributes
    """This class is a collective name space for basic calculations regarding Fanno flows.
        The constructor of this class can also determine the entire state of the flow given a
        partial state of the flow

    Args:
        gamma (float): ratio of specific heats
        mach (float, optional): mach number of the flow. Defaults to nan.
        t_t_st (float, optional): Ratio of Temperature to sonic temperature. Defaults to nan.
        p_p_st (float, optional): Ratio of Pressure to sonic pressure. Defaults to nan.
        rho_rho_st (float, optional): Ratio of density to sonic density. Defaults to nan.
        po_po_st (float, optional): Ratio of Total pressure to sonic total pressre. Defaults to nan.
        f4lst_d (float, optional): Effect of friction on pipe. Defaults to nan.
        u_u_st (float, optional): Velocity to sonic Velocity. Defaults to nan.
        flow_type (FlowState, optional):  States wether the flow is subsonic or supersonic. Used
                                        for Area Ratio Calculations. Defaults to FlowState.
                                        SUPER_SONIC.

    Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is
                                        given and flow state cannot be determined

    Useage:
        To use this class pass gamma and one of the known parameters of the flow and the
        rest are calculated.

    Valid_Combinations_of_Parameters:
        1: gamma, mach\n
        2: gamma, T/T*\n
        3: gamma, P/P*\n
        4: gamma, rho/rho*\n
        5: gamma, P0/P0*\n
        6: gamma, 4FL*/D, flowtype (flow type defaults to super sonic)\n
        7: gamma, U/U*\n
    """

    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        t_t_st: float = nan,
        p_p_st: float = nan,
        rho_rho_st: float = nan,
        po_po_st: float = nan,
        f4lst_d: float = nan,
        u_u_st: float = nan,
        flow_type: FlowState = FlowState.SUPER_SONIC,
    ) -> None:
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-statements

        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ Mach number of the flow """
        self.t_t_st = t_t_st
        """ Ratio of temperature to the sonic temperature T/T*"""
        self.p_p_st = p_p_st
        """ Ratio of pressure to the sonic pressure P/P*"""
        self.rho_rho_st = rho_rho_st
        """ Ratio of the density to the sonic density rho/rho* """
        self.po_po_st = po_po_st
        """ Ratio of the stagnation pressure to the sonic stagnation pressure P0/P0*"""
        self.f4lst_d = f4lst_d
        """ Friction parameter 4fL*/D"""
        self.u_u_st = u_u_st
        """ Ratio of velocity to the sonic velcoity U/U*"""
        self.flow_type = flow_type
        """ Type of flow which is either subsonic or supersonic (Type: flowstate) """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        # Pipe parameters
        self.choked_length = nan
        """ Pipe lenght required to choke the flow"""
        self.pipe_length = nan
        """ Current pipe lenght"""
        self.pipe_diameter = nan
        """ Diameter of the pipe"""
        self.friction_coeff = nan
        """ Friction coefficient"""

        # Down stream conditions if pipe parameters are given
        self.dwn_strm_mach = nan
        """ Mach number at the downstream point of the flow """
        self.dwn_strm_t_t_st = nan
        """ Ratio of temperature to the sonic temperature T/T* at the downstream point"""
        self.dwn_strm_p_p_st = nan
        """ Ratio of pressure to the sonic pressure P/P* at the downstream point"""
        self.dwn_strm_rho_rho_st = nan
        """ Ratio of the density to the sonic density rho/rho* at the downstream point"""
        self.dwn_strm_po_po_st = nan
        """ Ratio of stagnation pressure to sonic stagnation pressure P0/P0* at downstream point"""
        self.dwn_strm_f4lst_d = nan
        """ Friction parameter 4fL*/D of the downstream flow"""
        self.dwn_strm_u_u_st = nan
        """ Ratio of velocity to the sonic velcoity U/U* at the downstream point"""

        # Downstream / Initial Conditions
        self.t2_t1 = nan
        """ Ratio of downstream temperature to upstream temperature"""
        self.p2_p1 = nan
        """ Ratio of downstream pressure to upstream pressure"""
        self.rho2_rho1 = nan
        """ Ratio of downstream density to upstream density"""
        self.po2_po1 = nan
        """ Ratio of downstream total pressure to upstream total pressure"""
        self.f4ld2_f4ld1 = nan
        """ Ratio of the downstream friction parameter to the upstream friction parameter"""
        self.u2_u1 = nan
        """ Ratio of downstream velocity to upstream velocity """

        if not check_value(self.gamma):
            raise GammaNotDefinedError()

        if check_value(self.mach):
            pass

        elif check_value(self.t_t_st):
            self.mach = FannoFlowRelations.calc_mach_from_t_t_star(self.t_t_st, self.gamma)

        elif check_value(self.p_p_st):
            self.mach = FannoFlowRelations.calc_mach_from_p_p_star(self.p_p_st, self.gamma)

        elif check_value(self.rho_rho_st):
            self.mach = FannoFlowRelations.calc_mach_from_rho_rho_star(self.rho_rho_st, self.gamma)

        elif check_value(self.po_po_st):
            self.mach = FannoFlowRelations.calc_mach_from_po_po_star(
                self.po_po_st, self.gamma, flow_type=self.flow_type
            )

        elif check_value(self.f4lst_d, self.flow_type):
            self.mach = FannoFlowRelations.calc_mach_from_4flst_d(
                self.f4lst_d, self.gamma, flow_type=self.flow_type
            )

        elif check_value(self.u_u_st):
            self.mach = FannoFlowRelations.calc_mach_from_u_u_star(self.u_u_st, self.gamma)

        else:
            raise InvalidOptionCombinationError()

        if check_value(self.mach):
            self._calc_state()

    @property
    def choked_flow(self) -> bool:
        """True if the flow has reached the choked condition (I.E. Sonic)"""
        return bool(self.pipe_length >= self.choked_length)

    def __str__(self) -> str:
        color = Fore.GREEN if not self.choked_flow else Fore.YELLOW

        return "".join(
            [
                named_header("Fanno Relations at Mach", self.mach, precision=self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string("T/T*", self.t_t_st, self.precision, dot_line=True),
                to_string("P/P*", self.p_p_st, self.precision),
                to_string(f"{lcg.rho}/{lcg.rho}", self.rho_rho_st, self.precision, dot_line=True),
                to_string("4FL*/D", self.f4lst_d, self.precision),
                to_string("U/U*", self.u_u_st, self.precision, dot_line=True),
                to_string("Flow Type", self.flow_type.name, self.precision),
                seperator(),
                named_subheader("Pipe Parameters"),
                to_string("Length For Chocked Flow", self.choked_length, self.precision),
                to_string(
                    "Is Flow Choked? ", self.choked_flow, self.precision, dot_line=True, color=color
                ),
                to_string("Pipe Length", self.pipe_length, self.precision),
                to_string("Pipe Diameter", self.pipe_diameter, self.precision, dot_line=True),
                to_string("Friction Coefficient", self.friction_coeff, self.precision),
                seperator(),
                named_subheader("Down Stream Conditions"),
                to_string("Mach", self.dwn_strm_mach, self.precision),
                to_string("T/T*", self.t_t_st, self.precision, dot_line=True),
                to_string("P/P*", self.dwn_strm_p_p_st, self.precision),
                to_string("P0/P0*", self.dwn_strm_po_po_st, self.precision, dot_line=True),
                to_string(f"{lcg.rho}/{lcg.rho}", self.dwn_strm_rho_rho_st, self.precision),
                to_string("4FL*/D", self.dwn_strm_f4lst_d, self.precision, dot_line=True),
                to_string("U/U*", self.dwn_strm_u_u_st, self.precision),
                seperator(),
                named_subheader("Conditions Across Friction Area"),
                to_string("p2/p1", self.p2_p1, self.precision),
                to_string(
                    f"{lcg.rho}2/{lcg.rho}1",
                    self.rho2_rho1,
                    self.precision,
                    dot_line=True,
                ),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("4FL*/D2 / 4FL*/D 1", self.f4ld2_f4ld1, self.precision),
                to_string("U2/U1", self.u2_u1, self.precision, dot_line=True),
                footer(),
            ]
        )

    def apply_pipe_parameters(
        self, diameter: float, length: float, friction_coeff: float = 0.005
    ) -> None:
        """This functions applies parameters of a known pipe to the determined state of the flow.
            This allows the state at the downstream end of the pipe or pipe section to be found

        Args:
            diameter (float): Diameter of pipe
            length (float): Length of pipe
            friction_coeff (float, optional): Friction coefficient of the pipe. Defaults to 0.005
                                            which holds for Reynolds > 10^5.
        """
        self.pipe_diameter = diameter
        self.pipe_length = length
        self.friction_coeff = friction_coeff
        self._calc_down_stream_state()

    def _calc_state(self) -> None:
        self.t_t_st = FannoFlowRelations.calc_t_t_star(self.mach, self.gamma)
        self.p_p_st = FannoFlowRelations.calc_p_p_star(self.mach, self.gamma)
        self.rho_rho_st = FannoFlowRelations.calc_rho_rho_star(self.mach, self.gamma)
        self.po_po_st = FannoFlowRelations.calc_po_po_star(self.mach, self.gamma)
        self.f4lst_d = FannoFlowRelations.calc_4flst_d(self.mach, self.gamma)
        self.u_u_st = FannoFlowRelations.calc_u_u_starar(self.mach, self.gamma)

        if self.mach < 1:
            self.flow_type = FlowState.SUB_SONIC
        else:
            self.flow_type = FlowState.SUPER_SONIC

    def _calc_down_stream_state(self) -> None:
        if not check_value(self.pipe_diameter, self.pipe_length, self.friction_coeff):
            return

        self.choked_length = self.f4lst_d * self.pipe_diameter / 4 / self.friction_coeff
        self.dwn_strm_f4lst_d = (
            self.f4lst_d - 4 * self.friction_coeff * self.pipe_length / self.pipe_diameter
        )

        if self.pipe_length > self.choked_length:
            self.dwn_strm_mach = 1
        else:
            self.dwn_strm_mach = FannoFlowRelations.calc_mach_from_4flst_d(
                self.dwn_strm_f4lst_d, self.gamma, self.flow_type
            )

        self.dwn_strm_t_t_st = FannoFlowRelations.calc_t_t_star(self.dwn_strm_mach, self.gamma)
        self.dwn_strm_p_p_st = FannoFlowRelations.calc_p_p_star(self.dwn_strm_mach, self.gamma)
        self.dwn_strm_rho_rho_st = FannoFlowRelations.calc_rho_rho_star(
            self.dwn_strm_mach, self.gamma
        )
        self.dwn_strm_po_po_st = FannoFlowRelations.calc_po_po_star(self.dwn_strm_mach, self.gamma)
        self.dwn_strm_u_u_st = FannoFlowRelations.calc_u_u_starar(self.dwn_strm_mach, self.gamma)
        self.dwn_strm_f4lst_d = FannoFlowRelations.calc_4flst_d(self.dwn_strm_mach, self.gamma)

        # Calculate parameter ratios across the condition
        self.t2_t1 = self.dwn_strm_t_t_st / self.t_t_st
        self.p2_p1 = self.dwn_strm_p_p_st / self.p_p_st
        self.rho2_rho1 = self.dwn_strm_rho_rho_st / self.rho_rho_st
        self.po2_po1 = self.dwn_strm_po_po_st / self.po_po_st
        self.f4ld2_f4ld1 = self.dwn_strm_f4lst_d / self.f4lst_d
        self.u2_u1 = self.dwn_strm_u_u_st / self.u_u_st

    @staticmethod
    def calc_t_t_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static temperature to sonic temperature T/T*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: T/T*
        """
        return (gamma + 1) / (2 + (gamma - 1) * pow(mach, 2)) - offset

    @staticmethod
    def calc_mach_from_t_t_star(t_t_st: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of static temp to sonic static temp T/T*

        Args:
            t_t_st (float): Ratio of static temperature to sonic static temperature T/T*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """

        return brenth(
            FannoFlowRelations.calc_t_t_star,
            1e-9,
            40,
            args=(
                gamma,
                t_t_st,
            ),
        )  # type: ignore

    @staticmethod
    def calc_p_p_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static pressure to sonic pressure P/P*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: P/P*
        """
        return sqrt(FannoFlowRelations.calc_t_t_star(mach, gamma)) / mach - offset

    @staticmethod
    def calc_mach_from_p_p_star(p_p_st: float, gamma: float) -> float:
        """
        Calculates the mach based of the ratio of static pressure to sonic static pressure P/P*

        Args:
            p_p_st (float): Ratio of static pressure to sonic static pressure P/P*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(
            FannoFlowRelations.calc_p_p_star,
            1e-9,
            40,
            args=(
                gamma,
                p_p_st,
            ),
        )  # type: ignore

    @staticmethod
    def calc_rho_rho_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static density to sonic density Rho/Rho*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: Rho/Rho*
        """
        return sqrt(1 / FannoFlowRelations.calc_t_t_star(mach, gamma)) / mach - offset

    @staticmethod
    def calc_mach_from_rho_rho_star(rho_rho_st: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of density to sonic density Rho/Rho*

        Args:
            rho_rho_st (float): Ratio of density to sonic density Rho/Rho*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(
            FannoFlowRelations.calc_rho_rho_star,
            1e-9,
            40,
            args=(
                gamma,
                rho_rho_st,
            ),
        )  # type: ignore

    @staticmethod
    def calc_po_po_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static density to sonic density P0/P0*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: P0/P0*
        """
        gp1 = gamma + 1
        gm1 = gamma - 1
        return (
            pow(1 / FannoFlowRelations.calc_t_t_star(mach, gamma), gp1 / (2 * gm1)) / mach - offset
        )

    @staticmethod
    def calc_mach_from_po_po_star(
        po_po_st: float, gamma: float, flow_type: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach based of the ratio of total pressure to sonic total pressure P0/P0*

        Args:
            po_po_st (float): Ratio of total pressure to sonic total pressure P0/P0*
            gamma (float): ratio of specific heats
            flow_type (FlowState, optional): States whether the flow is currently supersonic or
                                            subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        tolerance = 1e-5
        if po_po_st == 1.0:
            return 1

        if flow_type == FlowState.SUPER_SONIC:
            return brenth(
                FannoFlowRelations.calc_po_po_star,
                1 + tolerance,
                40,
                args=(
                    gamma,
                    po_po_st,
                ),
            )  # type: ignore

        if flow_type == FlowState.SUB_SONIC:
            return brenth(
                FannoFlowRelations.calc_po_po_star,
                tolerance,
                1 - tolerance,
                args=(
                    gamma,
                    po_po_st,
                ),
            )  # type: ignore

        err = f"{Fore.RED} Flow Type [{flow_type}] not supported for Fanno Flow! {Fore.RESET}"
        raise ValueError(err)

    @staticmethod
    def calc_4flst_d(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates friction parameter for flow

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                        value. Defaults to 0.0.

        Returns:
            float: 4FL*/D
        """
        gp1 = gamma + 1
        t_t_st = FannoFlowRelations.calc_t_t_star(mach, gamma)
        m_sqr = pow(mach, 2)
        return (1 - m_sqr) / (gamma * m_sqr) + (gp1 / (2 * gamma)) * log(t_t_st * m_sqr) - offset

    @staticmethod
    def calc_mach_from_4flst_d(
        f4lst_d: float, gamma: float, flow_type: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach number from the friction parameter

        Args:
            f4lSt_d (float): friction parameter 4FL*/D
            gamma (float): ratio of specific heats
            flow_type (FlowState, optional):  Type of flow whether it is super sonic of subsonic.
                                            Defaults to FlowState.SUPER_SONIC.

        Raises:
            ValueError: Raised if Flow State is not supersonic of subsonic

        Returns:
            float: mach number
        """
        if f4lst_d == 0.0:
            return 1

        if flow_type == FlowState.SUPER_SONIC:
            return brenth(
                FannoFlowRelations.calc_4flst_d,
                1.00001,
                50,
                args=(
                    gamma,
                    f4lst_d,
                ),
            )  # type: ignore

        if flow_type == FlowState.SUB_SONIC:
            return brenth(
                FannoFlowRelations.calc_4flst_d,
                1e-5,
                0.9999,
                args=(
                    gamma,
                    f4lst_d,
                ),
            )  # type: ignore

        err = f"{Fore.RED}Flow Type [{flow_type}] not supported for FannoFlow!{Fore.RESET}"
        raise ValueError(err)

    @staticmethod
    def calc_u_u_starar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static velocity to sonic velocity U/U*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific
                                    value. Defaults to 0.0.

        Returns:
            float: U/U*
        """
        t_t_st = FannoFlowRelations.calc_t_t_star(mach, gamma)
        return mach * sqrt(t_t_st) - offset

    @staticmethod
    def calc_mach_from_u_u_star(u_u_st: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of velocity to sonic velocity U/U*

        Args:
            u_u_st (float): Ratio of velocity to sonic velocity U/U*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(
            FannoFlowRelations.calc_u_u_starar, 1e-9, 40, args=(gamma, u_u_st)
        )  # type: ignore
