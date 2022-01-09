from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore
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
        GammaNotDefinedError: [description]
        InvalidOptionCombinationError: [description]
        
    Useage:
        To use this class pass gamma and one of the known parameters of the flow and the rest are calculated. 

    Valid Combinations of Parameters:
        gamma, mach
        gamma, T/T*, flowtype (flow type defaults to super sonic)
        gamma, P/P*
        gamma, rho/rho*
        gamma, P0/P0*, flowtype (flow type defaults to super sonic)
        gamma, T0/T0*, flowtype (flow type defaults to super sonic)
        gamma, U/U*
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
        self.mach = mach
        self.t_tSt = t_tSt
        self.p_pSt = p_pSt
        self.rho_rhoSt = rho_rhoSt
        self.po_poSt = po_poSt
        self.u_uSt = u_uSt
        self.to_toSt = to_toSt
        self.flowType = flowType
        self.preciscion = 4

        # Pipe parameters
        self.chokedHeat = nan
        self.heat = nan
        self.gasConstantR = nan
        self.cp = nan

        # Down stream conditions if pipe parameters are given
        self.dwnStrmMach = nan
        self.dwnStrm_t_tSt = nan
        self.dwnStrm_p_pSt = nan
        self.dwnStrm_rho_rhoSt = nan
        self.dwnStrm_po_poSt = nan
        self.dwnStrm_u_uSt = nan
        self.dwnStrm_to_toSt = nan

        # Downstream / Initial Conditions
        self.t2_t1 = nan
        self.p2_p1 = nan
        self.rho2_rho1 = nan
        self.po2_po1 = nan
        self.to2_to1 = nan
        self.u2_u1 = nan
        self.to2 = nan
        self.to1 = nan

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
                named_header("Rayleigh Relations at Mach", self.mach, precision=self.preciscion),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.preciscion),
                to_string("T/T*", self.t_tSt, self.preciscion, dot_line=True),
                to_string("P/P*", self.p_pSt, self.preciscion),
                to_string("{}/{}*".format(lcg.rho, lcg.rho), self.rho_rhoSt, self.preciscion, dot_line=True),
                to_string("P0/P0*", self.po_poSt, self.preciscion),
                to_string("U/U*", self.u_uSt, self.preciscion, dot_line=True),
                to_string("T0/T0*", self.to_toSt, self.precision),
                to_string("Flow Type", self.flowType.name, self.preciscion, dot_line=True),
                seperator(),
                named_subheader("Pipe Parameters"),
                to_string("Heat Req. For Chocked Flow", self.chokedHeat, self.preciscion),
                color,
                to_string("Is Flow Choked? ", self.chockedFlow, self.preciscion, dot_line=True),
                to_string("Added Heat", self.heat, self.preciscion),
                to_string("Gas Constant R", self.gasConstantR, self.preciscion, dot_line=True),
                to_string("Cp", self.cp, self.preciscion),
                to_string("T01", self.to1, self.precision, dot_line=True),
                to_string("T02", self.to2, self.precision),
                seperator(),
                named_subheader("Down Stream Conditions"),
                to_string("Mach", self.dwnStrmMach, self.preciscion),
                to_string("T/T*", self.t_tSt, self.preciscion, dot_line=True),
                to_string("P/P*", self.dwnStrm_p_pSt, self.preciscion),
                to_string("P0/P0*", self.dwnStrm_po_poSt, self.preciscion, dot_line=True),
                to_string("{}/{}*".format(lcg.rho, lcg.rho), self.dwnStrm_rho_rhoSt, self.preciscion),
                to_string("T0/T0*", self.to_toSt, self.preciscion, dot_line=True),
                to_string("U/U*", self.dwnStrm_u_uSt, self.preciscion),
                seperator(),
                named_subheader("Conditions Across Heat Addition"),
                to_string("P2/P1", self.p2_p1, self.preciscion),
                to_string("{}2/{}1".format(lcg.rho, lcg.rho), self.rho2_rho1, self.preciscion, dot_line=True),
                to_string("T2/T1", self.t2_t1, self.preciscion),
                to_string("P02/P01", self.po2_po1, self.preciscion, dot_line=True),
                to_string("T02/T01", self.to2_to1, self.preciscion),
                to_string("U2/U1", self.u2_u1, self.preciscion, dot_line=True),
                footer(),
            ]
        )

    @staticmethod
    def calc_P_Pstar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates p/p* given gamma and mach"""
        return (1 + gamma) / (1 + gamma * pow(mach, 2)) - offset

    @staticmethod
    def calc_mach_from_P_PStar(p_pSt: float, gamma: float) -> float:
        """ Calculates mach given p/p* and gamma"""
        return brenth(RayleighFlowRelations.calc_P_Pstar, 1e-9, 40, args=(gamma, p_pSt,))

    @staticmethod
    def calc_T_Tstar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates T_T* given gamma and Mach, offset can be applied for root finding"""
        return pow(mach, 2) * pow(RayleighFlowRelations.calc_P_Pstar(mach, gamma), 2) - offset

    @staticmethod
    def calc_mach_from_T_TStar(
        t_tSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates Mach given a T_T* value and gamma"""
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
        """Calculates rho/rho* given gamma and mach"""
        return 1 / RayleighFlowRelations.calc_P_Pstar(mach, gamma) / pow(mach, 2) - offset

    @staticmethod
    def calc_mach_from_Rho_RhoStar(rho_rhoSt: float, gamma: float) -> float:
        """ Calculates mach given rho/rho* and gamma"""
        return brenth(RayleighFlowRelations.calc_Rho_RhoStar, 1e-9, 40, args=(gamma, rho_rhoSt,))

    @staticmethod
    def calc_Po_PoStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates po/po* given mach and gamma"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        p_pSt = RayleighFlowRelations.calc_P_Pstar(mach, gamma)
        ratio = (2 + gm1 * pow(mach, 2)) / gp1
        return p_pSt * pow(ratio, gamma / gm1) - offset

    @staticmethod
    def calc_mach_from_Po_PoStar(
        po_poSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """ Calculates mach given po/po* and gamma"""
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
        """Calculates Mach given a To_To* value and gamma"""
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
        """ Calculates U_U* given mach and gamma"""
        gp1 = gamma + 1
        mSqr = pow(mach, 2)
        return gp1 * mSqr / (1 + gamma * mSqr) - offset

    @staticmethod
    def calc_mach_from_U_USt(u_uSt: float, gamma: float) -> float:
        """ Calculates Mach given U_U* and gamma"""
        return brenth(RayleighFlowRelations.calc_U_UStar, 1e-9, 40, args=(gamma, u_uSt))
