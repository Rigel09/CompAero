from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from CompAero.internal import (
    FlowState,
    checkValue,
    value_to_string,
    named_subheader,
    named_header,
    footer,
    seperator,
)
from CompAero.greek_letters import LowerCaseGreek as lcg


class RayleighFlowRelations:
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
        self._preciscion = 4

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

        if checkValue(self.t_tSt):
            self.mach = RayleighFlowRelations.calcMachFrom_T_TSt(
                self.t_tSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.p_pSt):
            self.mach = RayleighFlowRelations.calcMachFrom_P_PSt(self.p_pSt, self.gamma)

        elif checkValue(self.rho_rhoSt):
            self.mach = RayleighFlowRelations.calcMachFrom_Rho_RhoSt(self.rho_rhoSt, self.gamma)

        elif checkValue(self.po_poSt):
            self.mach = RayleighFlowRelations.calcMachFrom_Po_PoSt(
                self.po_poSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.to_toSt):
            self.mach = RayleighFlowRelations.calcMachFrom_To_ToSt(
                self.to_toSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.u_uSt):
            self.mach = RayleighFlowRelations.calcMachFrom_U_USt(self.u_uSt, self.gamma)

        if checkValue(self.mach):
            self.__calculateState()

    @property
    def chockedFlow(self) -> bool:
        """ True if flow is choked """
        return self.heat > self.chokedHeat

    def __calculateState(self) -> None:
        self.t_tSt = RayleighFlowRelations.calcT_TSt_FromMach(self.mach, self.gamma)
        self.p_pSt = RayleighFlowRelations.calcP_PSt_FromMach(self.mach, self.gamma)
        self.rho_rhoSt = RayleighFlowRelations.calcRho_RhoSt_FromMach(self.mach, self.gamma)
        self.po_poSt = RayleighFlowRelations.calcPo_PoSt_FromMach(self.mach, self.gamma)
        self.to_toSt = RayleighFlowRelations.calcTo_ToSt_FromMach(self.mach, self.gamma)
        self.u_uSt = RayleighFlowRelations.calcU_USt_FromMach(self.mach, self.gamma)
        self.flowType = FlowState.SUPER_SONIC if self.mach > 1.0 else FlowState.SUB_SONIC

    def simulateHeatAddition(self, heat: float, to1: float, gasConstantR: float) -> None:
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

        self.dwnStrmMach = RayleighFlowRelations.calcMachFrom_To_ToSt(
            self.dwnStrm_to_toSt, self.gamma, self.flowType
        )
        self.dwnStrm_t_tSt = RayleighFlowRelations.calcT_TSt_FromMach(self.dwnStrmMach, self.gamma)
        self.dwnStrm_p_pSt = RayleighFlowRelations.calcP_PSt_FromMach(self.dwnStrmMach, self.gamma)
        self.dwnStrm_rho_rhoSt = RayleighFlowRelations.calcRho_RhoSt_FromMach(self.dwnStrmMach, self.gamma)
        self.dwnStrm_po_poSt = RayleighFlowRelations.calcPo_PoSt_FromMach(self.dwnStrmMach, self.gamma)
        self.dwnStrm_to_toSt = RayleighFlowRelations.calcTo_ToSt_FromMach(self.dwnStrmMach, self.gamma)
        self.dwnStrm_u_uSt = RayleighFlowRelations.calcU_USt_FromMach(self.dwnStrmMach, self.gamma)

        self.to2_to1 = self.dwnStrm_to_toSt / self.to_toSt
        self.t2_t1 = self.dwnStrm_t_tSt / self.t_tSt
        self.p2_p1 = self.dwnStrm_p_pSt / self.p_pSt
        self.rho2_rho1 = self.dwnStrm_rho_rhoSt / self.rho_rhoSt
        self.po2_po1 = self.dwnStrm_po_poSt / self.po_poSt
        self.u2_u1 = self.dwnStrm_u_uSt / self.u_uSt

    @property
    def precision(self) -> int:
        return self._preciscion

    @precision.setter
    def precision(self, precision) -> None:
        self._preciscion = precision

    def __str__(self) -> str:
        color = Back.GREEN + Fore.BLACK if self.heat < self.chokedHeat else Back.YELLOW + Fore.BLACK

        return "".join(
            [
                named_header("Rayleigh Relations at Mach", self.mach, precision=self._preciscion),
                seperator(),
                value_to_string(lcg.gamma, self.gamma, self._preciscion),
                value_to_string("T/T*", self.t_tSt, self._preciscion, dot_line=True),
                value_to_string("P/P*", self.p_pSt, self._preciscion),
                value_to_string(
                    "{}/{}*".format(lcg.rho, lcg.rho), self.rho_rhoSt, self._preciscion, dot_line=True
                ),
                value_to_string("P0/P0*", self.po_poSt, self._preciscion),
                value_to_string("U/U*", self.u_uSt, self._preciscion, dot_line=True),
                value_to_string("T0/T0*", self.to_toSt, self.precision),
                value_to_string("Flow Type", self.flowType.name, self._preciscion, dot_line=True),
                seperator(),
                named_subheader("Pipe Parameters"),
                value_to_string("Heat Req. For Chocked Flow", self.chokedHeat, self._preciscion),
                color,
                value_to_string("Is Flow Choked? ", self.chockedFlow, self._preciscion, dot_line=True),
                value_to_string("Added Heat", self.heat, self._preciscion),
                value_to_string("Gas Constant R", self.gasConstantR, self._preciscion, dot_line=True),
                value_to_string("Cp", self.cp, self._preciscion),
                value_to_string("T01", self.to1, self.precision, dot_line=True),
                value_to_string("T02", self.to2, self.precision),
                seperator(),
                named_subheader("Down Stream Conditions"),
                value_to_string("Mach", self.dwnStrmMach, self._preciscion),
                value_to_string("T/T*", self.t_tSt, self._preciscion, dot_line=True),
                value_to_string("P/P*", self.dwnStrm_p_pSt, self._preciscion),
                value_to_string("P0/P0*", self.dwnStrm_po_poSt, self._preciscion, dot_line=True),
                value_to_string("{}/{}*".format(lcg.rho, lcg.rho), self.dwnStrm_rho_rhoSt, self._preciscion),
                value_to_string("T0/T0*", self.to_toSt, self._preciscion, dot_line=True),
                value_to_string("U/U*", self.dwnStrm_u_uSt, self._preciscion),
                seperator(),
                named_subheader("Conditions Across Heat Addition"),
                value_to_string("P2/P1", self.p2_p1, self._preciscion),
                value_to_string(
                    "{}2/{}1".format(lcg.rho, lcg.rho), self.rho2_rho1, self._preciscion, dot_line=True
                ),
                value_to_string("T2/T1", self.t2_t1, self._preciscion),
                value_to_string("P02/P01", self.po2_po1, self._preciscion, dot_line=True),
                value_to_string("T02/T01", self.to2_to1, self._preciscion),
                value_to_string("U2/U1", self.u2_u1, self._preciscion, dot_line=True),
                footer(),
            ]
        )

    @staticmethod
    def calcP_PSt_FromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates p/p* given gamma and mach"""
        return (1 + gamma) / (1 + gamma * pow(mach, 2)) - offset

    @staticmethod
    def calcMachFrom_P_PSt(p_pSt: float, gamma: float) -> float:
        """ Calculates mach given p/p* and gamma"""
        return brenth(RayleighFlowRelations.calcP_PSt_FromMach, 1e-9, 40, args=(gamma, p_pSt,))

    @staticmethod
    def calcT_TSt_FromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates T_T* given gamma and Mach, offset can be applied for root finding"""
        return pow(mach, 2) * pow(RayleighFlowRelations.calcP_PSt_FromMach(mach, gamma), 2) - offset

    @staticmethod
    def calcMachFrom_T_TSt(t_tSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC) -> float:
        """Calculates Mach given a T_T* value and gamma"""
        tolerance = 1e-5
        if t_tSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(RayleighFlowRelations.calcT_TSt_FromMach, 1 + tolerance, 40, args=(gamma, t_tSt,))
        elif flowType == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calcT_TSt_FromMach, tolerance, 1 - tolerance, args=(gamma, t_tSt,)
            )
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calcRho_RhoSt_FromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates rho/rho* given gamma and mach"""
        return 1 / RayleighFlowRelations.calcP_PSt_FromMach(mach, gamma) / pow(mach, 2) - offset

    @staticmethod
    def calcMachFrom_Rho_RhoSt(rho_rhoSt: float, gamma: float) -> float:
        """ Calculates mach given rho/rho* and gamma"""
        return brenth(RayleighFlowRelations.calcRho_RhoSt_FromMach, 1e-9, 40, args=(gamma, rho_rhoSt,))

    @staticmethod
    def calcPo_PoSt_FromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates po/po* given mach and gamma"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        p_pSt = RayleighFlowRelations.calcP_PSt_FromMach(mach, gamma)
        ratio = (2 + gm1 * pow(mach, 2)) / gp1
        return p_pSt * pow(ratio, gamma / gm1) - offset

    @staticmethod
    def calcMachFrom_Po_PoSt(
        po_poSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """ Calculates mach given po/po* and gamma"""
        tolerance = 1e-5
        if po_poSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(
                RayleighFlowRelations.calcPo_PoSt_FromMach, 1 + tolerance, 40, args=(gamma, po_poSt,)
            )
        elif flowType == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calcPo_PoSt_FromMach, tolerance, 1 - tolerance, args=(gamma, po_poSt,)
            )
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calcTo_ToSt_FromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates To_To* given gamma and Mach, offset can be applied for root finding"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        mSqr = pow(mach, 2)
        return gp1 * mSqr / pow((1 + gamma * mSqr), 2) * (2 + gm1 * mSqr) - offset

    @staticmethod
    def calcMachFrom_To_ToSt(
        t_tSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates Mach given a To_To* value and gamma"""
        tolerance = 1e-5
        if t_tSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(RayleighFlowRelations.calcTo_ToSt_FromMach, 1 + tolerance, 40, args=(gamma, t_tSt,))
        elif flowType == FlowState.SUB_SONIC:
            return brenth(
                RayleighFlowRelations.calcTo_ToSt_FromMach, tolerance, 1 - tolerance, args=(gamma, t_tSt,)
            )
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calcU_USt_FromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates U_U* given mach and gamma"""
        gp1 = gamma + 1
        mSqr = pow(mach, 2)
        return gp1 * mSqr / (1 + gamma * mSqr) - offset

    @staticmethod
    def calcMachFrom_U_USt(u_uSt: float, gamma: float) -> float:
        """ Calculates Mach given U_U* and gamma"""
        return brenth(RayleighFlowRelations.calcU_USt_FromMach, 1e-9, 40, args=(gamma, u_uSt))
