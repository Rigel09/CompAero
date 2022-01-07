from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from CompAero.internal import FlowState, checkValue


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
        self.__preciscion = 4

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
        return self.__preciscion

    @precision.setter
    def precision(self, precision) -> None:
        self.__preciscion = precision

    def __str__(self) -> str:
        gammaStr = str(round(self.gamma, self.__preciscion))
        machStr = str(round(self.mach, self.__preciscion))
        tempStr = str(round(self.t_tSt, self.__preciscion))
        prStr = str(round(self.p_pSt, self.__preciscion))
        porStr = str(round(self.po_poSt, self.__preciscion))
        rhoStr = str(round(self.rho_rhoSt, self.__preciscion))
        UStr = str(round(self.u_uSt, self.__preciscion))
        toStr = str(round(self.to_toSt, self.__preciscion))

        # Pipe parameters
        chokedStr = str(round(self.chokedHeat, self.__preciscion))
        heatStr = str(round(self.heat, self.__preciscion))
        gasConstStr = str(round(self.gasConstantR, self.__preciscion))
        cpStr = str(round(self.cp, self.__preciscion))

        # Down stream conditions if pipe parameters are given
        dwnMachStr = str(round(self.dwnStrmMach, self.__preciscion))
        dwnTmpStr = str(round(self.dwnStrm_t_tSt, self.__preciscion))
        dwnPressStr = str(round(self.dwnStrm_p_pSt, self.__preciscion))
        dwnRhoStr = str(round(self.dwnStrm_rho_rhoSt, self.__preciscion))
        dwnPorStr = str(round(self.dwnStrm_po_poSt, self.__preciscion))
        dwnVelStr = str(round(self.dwnStrm_u_uSt, self.__preciscion))
        dwnTorStr = str(round(self.dwnStrm_to_toSt, self.__preciscion))

        # Downstream / Initial Conditions
        temp2Str = str(round(self.t2_t1, self.__preciscion))
        press2Str = str(round(self.p2_p1, self.__preciscion))
        rho2Str = str(round(self.rho2_rho1, self.__preciscion))
        po2Str = str(round(self.po2_po1, self.__preciscion))
        to2RStr = str(round(self.to2_to1, self.__preciscion))
        u2Str = str(round(self.u2_u1, self.__preciscion))
        to2Str = str(round(self.to2, self.__preciscion))
        to1Str = str(round(self.to1, self.__preciscion))

        width = 50 - 2  # Width minus edges
        sep = "|{:{width}}|\n".format("", width=width)

        # Initial Conditions
        header = "|{:=^{width}}|\n".format("Fanno Relations at Mach: {}".format(machStr), width=width)
        gamma = "|{:<{width}}{}|\n".format("\u03B3", gammaStr, width=width - len(gammaStr))
        tempR = "|{:-<{width}}{}|\n".format("T/T*", tempStr, width=width - len(tempStr))
        pressR = "|{:<{width}}{}|\n".format("P/P*", prStr, width=width - len(prStr))
        poR = "|{:-<{width}}{}|\n".format("Po/Po*", porStr, width=width - len(porStr))
        rhoR = "|{:<{width}}{}|\n".format("\u03C1/\u03C1*", rhoStr, width=width - len(rhoStr))
        uStR = "|{:-<{width}}{}|\n".format("U/U*", UStr, width=width - len(UStr))
        toR = "|{:<{width}}{}|\n".format("To/To*", toStr, width=width - len(toStr))
        flow = "|{:-<{width}}{}|\n".format("FlowType", self.flowType, width=width - len(self.flowType.name))

        # Pipe
        pipeSep = "|{:-^{width}}|\n".format("Pipe Conditions", width=width)
        if self.heat > self.chokedHeat:
            choked = (
                "|"
                + Fore.BLACK
                + Back.LIGHTYELLOW_EX
                + "{:<{width}}{}".format("Heat fro Choked Flow", chokedStr, width=width - len(chokedStr))
                + Style.RESET_ALL
                + "|\n"
            )
        else:
            choked = "|{:<{width}}{}|\n".format(
                "Heat fro Choked Flow", chokedStr, width=width - len(chokedStr)
            )
        heat = "|{:-<{width}}{}|\n".format("Added Heat", heatStr, width=width - len(heatStr))
        gasConst = "|{:<{width}}{}|\n".format("Gas Constant R", gasConstStr, width=width - len(gasConstStr))
        cp = "|{:-<{width}}{}|\n".format("Cp", cpStr, width=width - len(cpStr))
        to1 = "|{:<{width}}{}|\n".format("To1", to1Str, width=width - len(to1Str))
        to2 = "|{:-<{width}}{}|\n".format("To2", to2Str, width=width - len(to2Str))

        # Downstream
        dwnStrmSep = "|{:-^{width}}|\n".format("Down Stream Conditions", width=width)
        dwmMach = "|{:<{width}}{}|\n".format("Mach", dwnMachStr, width=width - len(dwnMachStr))
        dwnTmp = "|{:-<{width}}{}|\n".format("T/T*", dwnTmpStr, width=width - len(dwnTmpStr))
        dwnPress = "|{:<{width}}{}|\n".format("P/P*", dwnPressStr, width=width - len(dwnPressStr))
        dwnRho = "|{:-<{width}}{}|\n".format("\u03C1/\u03C1*", dwnRhoStr, width=width - len(dwnRhoStr))
        dwnPor = "|{:<{width}}{}|\n".format("Po/Po*", dwnPorStr, width=width - len(dwnPorStr))
        dwnVel = "|{:-<{width}}{}|\n".format("U/U*", dwnVelStr, width=width - len(dwnVelStr))
        dwnTor = "|{:<{width}}{}|\n".format("To/To*", dwnTorStr, width=width - len(dwnTorStr))

        # Jump Conditions
        jumpSep = "|{:-^{width}}|\n".format("Conditions Across Heat Addition", width=width)
        temp2 = "|{:<{width}}{}|\n".format("T2/T1", temp2Str, width=width - len(temp2Str))
        press2 = "|{:-<{width}}{}|\n".format("P2/P1", press2Str, width=width - len(press2Str))
        rho2 = "|{:<{width}}{}|\n".format("\u03C12/\u03C11", rho2Str, width=width - len(rho2Str))
        po2 = "|{:-<{width}}{}|\n".format("Po2/Po1", po2Str, width=width - len(po2Str))
        to2R = "|{:<{width}}{}|\n".format("To2/To1", to2RStr, width=width - len(to2RStr))
        u2 = "|{:-<{width}}{}|\n".format("U2/U1", u2Str, width=width - len(u2Str))

        finish = "|{:=^{width}}|\n".format("", width=width)

        return "".join(
            [
                header,
                sep,
                gamma,
                tempR,
                pressR,
                poR,
                rhoR,
                uStR,
                toR,
                flow,
                sep,
                pipeSep,
                choked,
                heat,
                gasConst,
                cp,
                to1,
                to2,
                sep,
                dwnStrmSep,
                dwmMach,
                dwnTmp,
                dwnPress,
                dwnRho,
                dwnPor,
                dwnVel,
                dwnTor,
                sep,
                jumpSep,
                temp2,
                press2,
                rho2,
                po2,
                to2R,
                u2,
                finish,
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
