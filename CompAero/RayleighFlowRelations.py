from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from enum import Enum
from CompAero.internal import (
    FlowState,
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    checkValue,
    to_string,
    named_subheader,
    named_header,
    footer,
    seperator,
)
from CompAero.greek_letters import LowerCaseGreek as lcg


RAYLEIGH_FLOW_VALID_OPTIONS = [
    "gamma, mach"
    "gamma, T/T*, flowtype"
    "gamma, P/P*"
    "gamma, rho/rho*"
    "gamma, P0/P0*, flowtype"
    "gamma, T0/T0*, flowtype"
    "gamma, U/U*"
]

class RayleighFlowChoice(Enum):
    GAMMA_MACH = "gamma, mach"
    GAMMA_T_T_ST = "gamma, T/T*, flowtype"
    GAMMA_P_P_ST = "gamma, P/P*"
    GAMMA_RHO_RHO_ST = "gamma, rho/rho*"
    GAMMA_PO_PO_ST = "gamma, P0/P0*, flowtype"
    GAMMA_4FLSTD_FLOW_TYPE = "gamma, T0/T0*, flowtype"
    GAMMA_U_U_ST = "gamma, U/U*"


class RayleighFlowRelations:
    """ This class is a collective name space for basic calculations regarding Fanno flows. 
        The constructor of this class can also determine the entire state of the flow given a partial state of the flow 

    Args:
        gamma (float): ratio of specific heats      
        mach (float, optional): mach number of the flow. Defaults to nan.
        t_tSt (float, optional): Ratio of Temperature to sonic temperature. Defaults to nan.
        p_pSt (float, optional): Ratio of Pressure to sonic pressure. Defaults to nan.
        rho_rhoSt (float, optional): Ratio of density to sonic density. Defaults to nan.
        po_poSt (float, optional): Ratio of Total pressure to sonic total pressre. Defaults to nan.
        to_toStar (float, optional): Ratio of Total Temperature to sonic total temperature. Defaults to nan.
        u_uSt (float, optional): Velocity to sonic Velocity. Defaults to nan.
        flowType (FlowState, optional):  States wether the flow is subsonic or supersonic. Used for Area Ratio Calculations. Defaults to FlowState.SUPER_SONIC.
        
    Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and flow state cannot be determined
        
    Useage:
        To use this class pass gamma and one of the known parameters of the flow and the rest are calculated. 

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
        t_tSt: float = nan,
        p_pSt: float = nan,
        rho_rhoSt: float = nan,
        po_poSt: float = nan,
        to_toSt: float = nan,
        u_uSt: float = nan,
        flowType: FlowState = FlowState.SUPER_SONIC,
    ) -> None:
        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ Mach number of the flow """
        self.t_tSt = t_tSt
        """ Ratio of temperature to the sonic temperature T/T*"""
        self.p_pSt = p_pSt
        """ Ratio of pressure to the sonic pressure P/P*"""
        self.rho_rhoSt = rho_rhoSt
        """ Ratio of the density to the sonic density rho/rho* """
        self.po_poSt = po_poSt
        """ Ratio of the stagnation pressure to the sonic stagnation pressure P0/P0*"""
        self.u_uSt = u_uSt
        """ Ratio of velocity to the sonic velcoity U/U*"""
        self.to_toSt = to_toSt
        """ Ratio of the total temperature to the sonic total temperature T0/T0*"""
        self.flowType = flowType
        """ Type of flow which is either subsonic or supersonic (Type: flowstate) """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        # Pipe parameters
        self.chokedHeat = nan
        """ Amount of heat required to choke the flow"""
        self.heat = nan
        """ Amount of heat added to the flow """
        self.gasConstantR = nan
        """ Gas constant used for flow calculations"""
        self.cp = nan
        """ Specific heat Coefficient at constant pressure for the gas used"""

        # Down stream conditions if pipe parameters are given
        self.dwnStrmMach = nan
        """ Mach number at the downstream point of the flow """
        self.dwnStrm_t_tSt = nan
        """ Ratio of temperature to the sonic temperature T/T* at the downstream point"""
        self.dwnStrm_p_pSt = nan
        """ Ratio of pressure to the sonic pressure P/P* at the downstream point"""
        self.dwnStrm_rho_rhoSt = nan
        """ Ratio of the density to the sonic density rho/rho* at the downstream point"""
        self.dwnStrm_po_poSt = nan
        """ Ratio of the stagnation pressure to the sonic stagnation pressure P0/P0* at the downstream point"""
        self.dwnStrm_u_uSt = nan
        """ Ratio of velocity to the sonic velcoity U/U* at the downstream point"""
        self.dwnStrm_to_toSt = nan
        """ Ratio of the total temperature to the sonic total temperature T0/T0* at the downstream point"""

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

        if not checkValue(self.gamma):
            raise GammaNotDefinedError()

        if checkValue(self.mach):
            pass

        elif checkValue(self.t_tSt):
            self.mach = RayleighFlowRelations.calc_mach_from_T_TStar(
                self.t_tSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.p_pSt):
            self.mach = RayleighFlowRelations.calc_mach_from_P_PStar(self.p_pSt, self.gamma)

        elif checkValue(self.rho_rhoSt):
            self.mach = RayleighFlowRelations.calc_mach_from_Rho_RhoStar(self.rho_rhoSt, self.gamma)

        elif checkValue(self.po_poSt):
            self.mach = RayleighFlowRelations.calc_mach_from_Po_PoStar(
                self.po_poSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.to_toSt):
            self.mach = RayleighFlowRelations.calc_mach_from_To_ToStar(
                self.to_toSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.u_uSt):
            self.mach = RayleighFlowRelations.calc_mach_from_U_USt(self.u_uSt, self.gamma)

        else:
            raise InvalidOptionCombinationError()

        if checkValue(self.mach):
            self.__calculateState()

    @property
    def chockedFlow(self) -> bool:
        """ True if the added heat is greater than the heat addition required to choke the flow """
        return self.heat > self.chokedHeat

    def __calculateState(self) -> None:
        self.t_tSt = RayleighFlowRelations.calc_T_Tstar(self.mach, self.gamma)
        self.p_pSt = RayleighFlowRelations.calc_P_Pstar(self.mach, self.gamma)
        self.rho_rhoSt = RayleighFlowRelations.calc_Rho_RhoStar(self.mach, self.gamma)
        self.po_poSt = RayleighFlowRelations.calc_Po_PoStar(self.mach, self.gamma)
        self.to_toSt = RayleighFlowRelations.calc_To_ToSt(self.mach, self.gamma)
        self.u_uSt = RayleighFlowRelations.calc_U_UStar(self.mach, self.gamma)
        self.flowType = FlowState.SUPER_SONIC if self.mach > 1.0 else FlowState.SUB_SONIC

    def simulate_heat_addition(self, heat: float, to1: float, gasConstantR: float) -> None:
        self.gasConstantR = gasConstantR
        self.to1 = to1
        self.heat = heat
        self.cp = self.gamma * gasConstantR / (self.gamma - 1)
        self.to2 = self.to1 + self.heat / self.cp

        tempStar = self.to1 / self.to_toSt
        self.chokedHeat = self.cp * (tempStar - self.to1)

        self.to2_to1 = self.to2 / self.to1
        self.dwnStrm_to_toSt = self.to2_to1 * self.to_toSt

        if self.heat > self.chokedHeat:
            self.dwnStrm_to_toSt = 1

        self.dwnStrmMach = RayleighFlowRelations.calc_mach_from_To_ToStar(
            self.dwnStrm_to_toSt, self.gamma, self.flowType
        )
        self.dwnStrm_t_tSt = RayleighFlowRelations.calc_T_Tstar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_p_pSt = RayleighFlowRelations.calc_P_Pstar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_rho_rhoSt = RayleighFlowRelations.calc_Rho_RhoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_po_poSt = RayleighFlowRelations.calc_Po_PoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_to_toSt = RayleighFlowRelations.calc_To_ToSt(self.dwnStrmMach, self.gamma)
        self.dwnStrm_u_uSt = RayleighFlowRelations.calc_U_UStar(self.dwnStrmMach, self.gamma)

        self.to2_to1 = self.dwnStrm_to_toSt / self.to_toSt
        self.t2_t1 = self.dwnStrm_t_tSt / self.t_tSt
        self.p2_p1 = self.dwnStrm_p_pSt / self.p_pSt
        self.rho2_rho1 = self.dwnStrm_rho_rhoSt / self.rho_rhoSt
        self.po2_po1 = self.dwnStrm_po_poSt / self.po_poSt
        self.u2_u1 = self.dwnStrm_u_uSt / self.u_uSt

    def __str__(self) -> str:
        color = Back.GREEN + Fore.BLACK if self.heat < self.chokedHeat else Back.YELLOW + Fore.BLACK

        return "".join(
            [
                named_header("Rayleigh Relations at Mach", self.mach, precision=self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string("T/T*", self.t_tSt, self.precision, dot_line=True),
                to_string("P/P*", self.p_pSt, self.precision),
                to_string("{}/{}*".format(lcg.rho, lcg.rho), self.rho_rhoSt, self.precision, dot_line=True),
                to_string("P0/P0*", self.po_poSt, self.precision),
                to_string("U/U*", self.u_uSt, self.precision, dot_line=True),
                to_string("T0/T0*", self.to_toSt, self.precision),
                to_string("Flow Type", self.flowType.name, self.precision, dot_line=True),
                seperator(),
                named_subheader("Pipe Parameters"),
                to_string("Heat Req. For Chocked Flow", self.chokedHeat, self.precision),
                color,
                to_string("Is Flow Choked? ", self.chockedFlow, self.precision, dot_line=True),
                to_string("Added Heat", self.heat, self.precision),
                to_string("Gas Constant R", self.gasConstantR, self.precision, dot_line=True),
                to_string("Cp", self.cp, self.precision),
                to_string("T01", self.to1, self.precision, dot_line=True),
                to_string("T02", self.to2, self.precision),
                seperator(),
                named_subheader("Down Stream Conditions"),
                to_string("Mach", self.dwnStrmMach, self.precision),
                to_string("T/T*", self.t_tSt, self.precision, dot_line=True),
                to_string("P/P*", self.dwnStrm_p_pSt, self.precision),
                to_string("P0/P0*", self.dwnStrm_po_poSt, self.precision, dot_line=True),
                to_string("{}/{}*".format(lcg.rho, lcg.rho), self.dwnStrm_rho_rhoSt, self.precision),
                to_string("T0/T0*", self.to_toSt, self.precision, dot_line=True),
                to_string("U/U*", self.dwnStrm_u_uSt, self.precision),
                seperator(),
                named_subheader("Conditions Across Heat Addition"),
                to_string("P2/P1", self.p2_p1, self.precision),
                to_string("{}2/{}1".format(lcg.rho, lcg.rho), self.rho2_rho1, self.precision, dot_line=True),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("T02/T01", self.to2_to1, self.precision),
                to_string("U2/U1", self.u2_u1, self.precision, dot_line=True),
                footer(),
            ]
        )

    @staticmethod
    def calc_P_Pstar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static pressure to sonic pressure P/P*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: P/P*
        """
        return (1 + gamma) / (1 + gamma * pow(mach, 2)) - offset

    @staticmethod
    def calc_mach_from_P_PStar(p_pSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of static pressure to sonic static pressure P/P*

        Args:
            p_pSt (float): Ratio of static pressure to sonic static pressure P/P*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(RayleighFlowRelations.calc_P_Pstar, 1e-9, 40, args=(gamma, p_pSt,))

    @staticmethod
    def calc_T_Tstar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static temperature to sonic temperature T/T*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: T/T*
        """
        return pow(mach, 2) * pow(RayleighFlowRelations.calc_P_Pstar(mach, gamma), 2) - offset

    @staticmethod
    def calc_mach_from_T_TStar(
        t_tSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach number based of the ratio of static temperature to sonic static temperature T/T*

        Args:
            t_tSt (float): Ratio of static temperature to sonic static temperature T/T*
            gamma (float): ratio of specific heats
            flowType (FlowState, optional): States whether the flow is currently supersonic or subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        tolerance = 1e-5
        if t_tSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(RayleighFlowRelations.calc_T_Tstar, 1 + tolerance, 40, args=(gamma, t_tSt,))
        elif flowType == FlowState.SUB_SONIC:
            return brenth(RayleighFlowRelations.calc_T_Tstar, tolerance, 1 - tolerance, args=(gamma, t_tSt,))
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calc_Rho_RhoStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static density to sonic density Rho/Rho*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: Rho/Rho*
        """
        return 1 / RayleighFlowRelations.calc_P_Pstar(mach, gamma) / pow(mach, 2) - offset

    @staticmethod
    def calc_mach_from_Rho_RhoStar(rho_rhoSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of density to sonic density Rho/Rho*

        Args:
            rho_rhoSt (float): Ratio of density to sonic density Rho/Rho*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(RayleighFlowRelations.calc_Rho_RhoStar, 1e-9, 40, args=(gamma, rho_rhoSt,))

    @staticmethod
    def calc_Po_PoStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static density to sonic density P0/P0*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: P0/P0*
        """
        gp1 = gamma + 1
        gm1 = gamma - 1
        p_pSt = RayleighFlowRelations.calc_P_Pstar(mach, gamma)
        ratio = (2 + gm1 * pow(mach, 2)) / gp1
        return p_pSt * pow(ratio, gamma / gm1) - offset

    @staticmethod
    def calc_mach_from_Po_PoStar(
        po_poSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach number based of the ratio of total pressure to sonic total pressure P0/P0*

        Args:
            po_poSt (float): Ratio of total pressure to sonic total pressure P0/P0*
            gamma (float): ratio of specific heats
            flowType (FlowState, optional): States whether the flow is currently supersonic or subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        tolerance = 1e-5
        if po_poSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(RayleighFlowRelations.calc_Po_PoStar, 1 + tolerance, 40, args=(gamma, po_poSt,))
        elif flowType == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calc_Po_PoStar, tolerance, 1 - tolerance, args=(gamma, po_poSt,)
            )
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calc_To_ToSt(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates To_To* given gamma and Mach, offset can be applied for root finding"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        mSqr = pow(mach, 2)
        return gp1 * mSqr / pow((1 + gamma * mSqr), 2) * (2 + gm1 * mSqr) - offset

    @staticmethod
    def calc_mach_from_To_ToStar(
        t_tSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach number based of the ratio of total temperature to sonic total temperature T0/T0*

        Args:
            po_poSt (float): Ratio of total temperature to sonic total temperature T0/T0*
            gamma (float): ratio of specific heats
            flowType (FlowState, optional): States whether the flow is currently supersonic or subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        tolerance = 1e-5
        if t_tSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(RayleighFlowRelations.calc_To_ToSt, 1 + tolerance, 40, args=(gamma, t_tSt,))
        elif flowType == FlowState.SUB_SONIC:
            return brenth(RayleighFlowRelations.calc_To_ToSt, tolerance, 1 - tolerance, args=(gamma, t_tSt,))
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calc_U_UStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of static velocity to sonic velocity U/U*

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: U/U*
        """
        gp1 = gamma + 1
        mSqr = pow(mach, 2)
        return gp1 * mSqr / (1 + gamma * mSqr) - offset

    @staticmethod
    def calc_mach_from_U_USt(u_uSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of velocity to sonic velocity U/U*

        Args:
            u_uSt (float): Ratio of velocity to sonic velocity U/U*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(RayleighFlowRelations.calc_U_UStar, 1e-9, 40, args=(gamma, u_uSt))
