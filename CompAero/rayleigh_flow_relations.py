"""
This module contains the relationships for rayleigh flows
"""

from enum import Enum
from math import nan

from colorama import Back, Fore
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


class RayleighFlowChoice(Enum):
    """Valid choices for Rayleigh Flow initialization"""

    GAMMA_MACH = "gamma, mach"
    GAMMA_T_T_ST = "gamma, T/T*, flowtype"
    GAMMA_P_P_ST = "gamma, P/P*"
    GAMMA_RHO_RHO_ST = "gamma, rho/rho*"
    GAMMA_PO_PO_ST = "gamma, P0/P0*, flowtype"
    GAMMA_TO_TO_FLOW_TYPE = "gamma, T0/T0*, flowtype"
    GAMMA_U_U_ST = "gamma, U/U*"


RAYLEIGH_FLOW_VALID_OPTIONS = [x.value for x in RayleighFlowChoice]


class RayleighFlowRelations:
    # pylint: disable=too-many-instance-attributes
    """
    This class is a collective name space for basic calculations regarding Fanno flows.

    Args:
        gamma (float): ratio of specific heats
        mach (float, optional): mach number of the flow. Defaults to nan.
        t_t_st (float, optional): Ratio of Temperature to sonic temperature. Defaults to nan.
        p_p_st (float, optional): Ratio of Pressure to sonic pressure. Defaults to nan.
        rho_rho_st (float, optional): Ratio of density to sonic density. Defaults to nan.
        po_po_st (float, optional): Total pressure to sonic total pressure ratio. Defaults to nan.
        to_to_star (float, optional): Ratio of Total Tempto sonic total temp. Defaults to nan.
        u_u_st (float, optional): Velocity to sonic Velocity. Defaults to nan.
        flow_type (FlowState, optional):  States wether the flow is subsonic or supersonic. Used
                                        for Area Ratio Calculations. Defaults to
                                        FlowState.SUPER_SONIC.

    Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given
                                        and flow state cannot be determined

    Useage:
        To use this class pass gamma and one of the known parameters of the flow

    Valid_Combinations_of_Parameters:
        1: gamma, mach\n
        2: gamma, T/T*, flowtype (flow type defaults to super sonic)\n
        3: gamma, P/P*\n
        4: gamma, rho/rho*\n
        5: gamma, P0/P0*, flowtype (flow type defaults to super sonic)\n
        6: gamma, T0/T0*, flowtype (flow type defaults to super sonic)\n
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
        to_to_st: float = nan,
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
        self.u_u_st = u_u_st
        """ Ratio of velocity to the sonic velcoity U/U*"""
        self.to_to_st = to_to_st
        """ Ratio of the total temperature to the sonic total temperature T0/T0*"""
        self.flow_type = flow_type
        """ Type of flow which is either subsonic or supersonic (Type: flowstate) """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        # Pipe parameters
        self.choked_heat = nan
        """ Amount of heat required to choke the flow"""
        self.heat = nan
        """ Amount of heat added to the flow """
        self.gas_constant_r = nan
        """ Gas constant used for flow calculations"""
        self.cp = nan
        """ Specific heat Coefficient at constant pressure for the gas used"""

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
        """ Stagnation to sonic stagnation pressure ratio P0/P0* at the down stream point"""
        self.dwn_strm_u_u_st = nan
        """ Ratio of velocity to the sonic velcoity U/U* at the downstream point"""
        self.dwn_strm_to_to_st = nan
        """ Ratio of the total temp to the sonic total temp T0/T0* at the downstream point"""

        # Downstream / Initial Conditions
        self.t2_t1 = nan
        """ Ratio of downstream temperature to upstream temperature"""
        self.p2_p1 = nan
        """ Ratio of downstream pressure to upstream pressure"""
        self.rho2_rho1 = nan
        """ Ratio of downstream density to upstream density"""
        self.po2_po1 = nan
        """ Ratio of downstream total pressure to upstream total pressure"""
        self.to2_to1 = nan
        """ Ratio of downstream total temperature to upstream total temperature"""
        self.u2_u1 = nan
        """ Ratio of downstream velocity to upstream velocity """
        self.to2 = nan
        """ Downstream total temperature """
        self.to1 = nan
        """ Upstream total temperature"""

        if not check_value(self.gamma):
            raise GammaNotDefinedError()

        if check_value(self.mach):
            pass

        elif check_value(self.t_t_st):
            self.mach = RayleighFlowRelations.calc_mach_from_t_t_star(
                self.t_t_st, self.gamma, flow_type=self.flow_type
            )

        elif check_value(self.p_p_st):
            self.mach = RayleighFlowRelations.calc_mach_from_p_p_star(self.p_p_st, self.gamma)

        elif check_value(self.rho_rho_st):
            self.mach = RayleighFlowRelations.calc_mach_from_rho_rho_star(
                self.rho_rho_st, self.gamma
            )

        elif check_value(self.po_po_st):
            self.mach = RayleighFlowRelations.calc_mach_from_po_po_star(
                self.po_po_st, self.gamma, flow_type=self.flow_type
            )

        elif check_value(self.to_to_st):
            self.mach = RayleighFlowRelations.calc_mach_from_to_to_star(
                self.to_to_st, self.gamma, flow_type=self.flow_type
            )

        elif check_value(self.u_u_st):
            self.mach = RayleighFlowRelations.calc_mach_from_u_u_star(self.u_u_st, self.gamma)

        else:
            raise InvalidOptionCombinationError()

        if check_value(self.mach):
            self._calc_state()

    @property
    def choked_flow(self) -> bool:
        """True if the added heat is greater than the heat addition required to choke the flow"""
        return self.heat > self.choked_heat

    def _calc_state(self) -> None:
        self.t_t_st = RayleighFlowRelations.calc_t_t_star(self.mach, self.gamma)
        self.p_p_st = RayleighFlowRelations.calc_p_p_star(self.mach, self.gamma)
        self.rho_rho_st = RayleighFlowRelations.calc_rho_rho_star(self.mach, self.gamma)
        self.po_po_st = RayleighFlowRelations.calc_po_po_star(self.mach, self.gamma)
        self.to_to_st = RayleighFlowRelations.calc_to_to_star(self.mach, self.gamma)
        self.u_u_st = RayleighFlowRelations.calc_u_u_starar(self.mach, self.gamma)
        self.flow_type = FlowState.SUPER_SONIC if self.mach > 1.0 else FlowState.SUB_SONIC

    def simulate_heat_addition(self, heat: float, to1: float, gas_const_r: float) -> None:
        """
        This functions calculates the downstream flow conditions caused due to applied heat.

        Args:
            heat (float): The amount of heat being applied
            to1 (float): The initial temperature of the flow
            gas_const_r (float): The gas constant of the flow
        """
        self.gas_constant_r = gas_const_r
        self.to1 = to1
        self.heat = heat
        self.cp = self.gamma * gas_const_r / (self.gamma - 1)
        self.to2 = self.to1 + self.heat / self.cp

        temp_star = self.to1 / self.to_to_st
        self.choked_heat = self.cp * (temp_star - self.to1)

        self.to2_to1 = self.to2 / self.to1
        self.dwn_strm_to_to_st = self.to2_to1 * self.to_to_st

        if self.heat > self.choked_heat:
            self.dwn_strm_to_to_st = 1

        self.dwn_strm_mach = RayleighFlowRelations.calc_mach_from_to_to_star(
            self.dwn_strm_to_to_st, self.gamma, self.flow_type
        )
        self.dwn_strm_t_t_st = RayleighFlowRelations.calc_t_t_star(self.dwn_strm_mach, self.gamma)
        self.dwn_strm_p_p_st = RayleighFlowRelations.calc_p_p_star(self.dwn_strm_mach, self.gamma)
        self.dwn_strm_rho_rho_st = RayleighFlowRelations.calc_rho_rho_star(
            self.dwn_strm_mach, self.gamma
        )
        self.dwn_strm_po_po_st = RayleighFlowRelations.calc_po_po_star(
            self.dwn_strm_mach, self.gamma
        )
        self.dwn_strm_to_to_st = RayleighFlowRelations.calc_to_to_star(
            self.dwn_strm_mach, self.gamma
        )
        self.dwn_strm_u_u_st = RayleighFlowRelations.calc_u_u_starar(self.dwn_strm_mach, self.gamma)

        self.to2_to1 = self.dwn_strm_to_to_st / self.to_to_st
        self.t2_t1 = self.dwn_strm_t_t_st / self.t_t_st
        self.p2_p1 = self.dwn_strm_p_p_st / self.p_p_st
        self.rho2_rho1 = self.dwn_strm_rho_rho_st / self.rho_rho_st
        self.po2_po1 = self.dwn_strm_po_po_st / self.po_po_st
        self.u2_u1 = self.dwn_strm_u_u_st / self.u_u_st

    def __str__(self) -> str:
        color = (
            Back.GREEN + Fore.BLACK if self.heat < self.choked_heat else Back.YELLOW + Fore.BLACK
        )

        return "".join(
            [
                named_header("Rayleigh Relations at Mach", self.mach, precision=self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string("T/T*", self.t_t_st, self.precision, dot_line=True),
                to_string("P/P*", self.p_p_st, self.precision),
                to_string(f"{lcg.rho}/{lcg.rho}", self.rho_rho_st, self.precision, dot_line=True),
                to_string("P0/P0*", self.po_po_st, self.precision),
                to_string("U/U*", self.u_u_st, self.precision, dot_line=True),
                to_string("T0/T0*", self.to_to_st, self.precision),
                to_string("Flow Type", self.flow_type.name, self.precision, dot_line=True),
                seperator(),
                named_subheader("Pipe Parameters"),
                to_string("Heat Req. For Chocked Flow", self.choked_heat, self.precision),
                color,
                to_string("Is Flow Choked? ", self.choked_flow, self.precision, dot_line=True),
                to_string("Added Heat", self.heat, self.precision),
                to_string("Gas Constant R", self.gas_constant_r, self.precision, dot_line=True),
                to_string("Cp", self.cp, self.precision),
                to_string("T01", self.to1, self.precision, dot_line=True),
                to_string("T02", self.to2, self.precision),
                seperator(),
                named_subheader("Down Stream Conditions"),
                to_string("Mach", self.dwn_strm_mach, self.precision),
                to_string("T/T*", self.t_t_st, self.precision, dot_line=True),
                to_string("P/P*", self.dwn_strm_p_p_st, self.precision),
                to_string("P0/P0*", self.dwn_strm_po_po_st, self.precision, dot_line=True),
                to_string(f"{lcg.rho}/{lcg.rho}*", self.dwn_strm_rho_rho_st, self.precision),
                to_string("T0/T0*", self.to_to_st, self.precision, dot_line=True),
                to_string("U/U*", self.dwn_strm_u_u_st, self.precision),
                seperator(),
                named_subheader("Conditions Across Heat Addition"),
                to_string("P2/P1", self.p2_p1, self.precision),
                to_string(
                    f"{lcg.rho}2/{lcg.rho}1",
                    self.rho2_rho1,
                    self.precision,
                    dot_line=True,
                ),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("T02/T01", self.to2_to1, self.precision),
                to_string("U2/U1", self.u2_u1, self.precision, dot_line=True),
                footer(),
            ]
        )

    @staticmethod
    def calc_p_p_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static pressure to sonic pressure P/P*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding. Defaults to 0.0.

        Returns:
            float: P/P*
        """
        return (1 + gamma) / (1 + gamma * pow(mach, 2)) - offset

    @staticmethod
    def calc_mach_from_p_p_star(p_p_st: float, gamma: float) -> float:
        """
        Calculates the mach number based of the ratio of static pressure to sonic static
        pressure P/P*

        Args:
            p_p_st (float): Ratio of static pressure to sonic static pressure P/P*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(
            RayleighFlowRelations.calc_p_p_star,
            1e-9,
            40,
            args=(
                gamma,
                p_p_st,
            ),
        )  # type: ignore

    @staticmethod
    def calc_t_t_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static temperature to sonic temperature T/T*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding. Defaults to 0.0.

        Returns:
            float: T/T*
        """
        return pow(mach, 2) * pow(RayleighFlowRelations.calc_p_p_star(mach, gamma), 2) - offset

    @staticmethod
    def calc_mach_from_t_t_star(
        t_t_st: float, gamma: float, flow_type: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """
        Calculates the mach number based of the ratio of static tempto sonic static temp T/T*

        Args:
            t_t_st (float): Ratio of static temperature to sonic static temperature T/T*
            gamma (float): ratio of specific heats
            flow_type (FlowState, optional): States whether the flow is currently supersonic or
                                            subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        tolerance = 1e-5
        if t_t_st == 1.0:
            return 1

        if flow_type == FlowState.SUPER_SONIC:
            return brenth(
                RayleighFlowRelations.calc_t_t_star,
                1 + tolerance,
                40,
                args=(
                    gamma,
                    t_t_st,
                ),
            )  # type: ignore

        if flow_type == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calc_t_t_star,
                tolerance,
                1 - tolerance,
                args=(
                    gamma,
                    t_t_st,
                ),
            )  # type: ignore

        raise ValueError(
            f"{Fore.RED} Flow Type [{flow_type}] not supported for Fanno Flow!{Fore.RESET}"
        )

    @staticmethod
    def calc_rho_rho_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static density to sonic density Rho/Rho*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding. Defaults to 0.0.

        Returns:
            float: Rho/Rho*
        """
        return 1 / RayleighFlowRelations.calc_p_p_star(mach, gamma) / pow(mach, 2) - offset

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
            RayleighFlowRelations.calc_rho_rho_star,
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
            offset (float, optional): offset that can be used for root finding. Defaults to 0.0.

        Returns:
            float: P0/P0*
        """
        gp1 = gamma + 1
        gm1 = gamma - 1
        p_p_st = RayleighFlowRelations.calc_p_p_star(mach, gamma)
        ratio = (2 + gm1 * pow(mach, 2)) / gp1
        return p_p_st * pow(ratio, gamma / gm1) - offset

    @staticmethod
    def calc_mach_from_po_po_star(
        po_po_st: float, gamma: float, flow_type: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """
        Calculates the mach number based of the ratio of total pressure to sonic total pressure
        P0/P0*

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
                RayleighFlowRelations.calc_po_po_star,
                1 + tolerance,
                40,
                args=(
                    gamma,
                    po_po_st,
                ),
            )  # type: ignore

        if flow_type == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calc_po_po_star,
                tolerance,
                1 - tolerance,
                args=(
                    gamma,
                    po_po_st,
                ),
            )  # type: ignore

        raise ValueError(
            f"{Fore.RED} Flow Type [{flow_type}] not supported for Fanno Flow!{Fore.RESET}"
        )

    @staticmethod
    def calc_to_to_star(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates To_To* given gamma and Mach, offset can be applied for root finding"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        m_sqr = pow(mach, 2)
        return gp1 * m_sqr / pow((1 + gamma * m_sqr), 2) * (2 + gm1 * m_sqr) - offset

    @staticmethod
    def calc_mach_from_to_to_star(
        t_t_st: float, gamma: float, flow_type: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """
        Calculates the mach number based of the ratio of total tempto sonic total temp T0/T0*

        Args:
            po_po_st (float): Ratio of total temperature to sonic total temperature T0/T0*
            gamma (float): ratio of specific heats
            flow_type (FlowState, optional): States whether the flow is currently supersonic or
                                                subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        tolerance = 1e-5
        if t_t_st == 1.0:
            return 1

        if flow_type == FlowState.SUPER_SONIC:
            return brenth(
                RayleighFlowRelations.calc_to_to_star,
                1 + tolerance,
                40,
                args=(
                    gamma,
                    t_t_st,
                ),
            )  # type: ignore

        if flow_type == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calc_to_to_star,
                tolerance,
                1 - tolerance,
                args=(
                    gamma,
                    t_t_st,
                ),
            )  # type: ignore

        raise ValueError(
            f"{Fore.RED} Flow Type [{flow_type}] not supported for Fanno Flow!{Fore.RESET}"
        )

    @staticmethod
    def calc_u_u_starar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static velocity to sonic velocity U/U*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding. Defaults to 0.0.

        Returns:
            float: U/U*
        """
        gp1 = gamma + 1
        m_sqr = pow(mach, 2)
        return gp1 * m_sqr / (1 + gamma * m_sqr) - offset

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
            RayleighFlowRelations.calc_u_u_starar, 1e-9, 40, args=(gamma, u_u_st)
        )  # type: ignore
