from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore


class FannoFlowRelations:
    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        t_tSt: float = nan,
        p_pSt: float = nan,
        rho_rhoSt: float = nan,
        po_poSt: float = nan,
        f4LSt_D: float = nan,
        u_uSt: float = nan,
        flowType: str = "Supersonic",
    ) -> None:
        self.gamma = gamma
        self.mach = mach
        self.t_tSt = t_tSt
        self.p_pSt = p_pSt
        self.rho_rhoSt = rho_rhoSt
        self.po_poSt = po_poSt
        self.f4LSt_D = f4LSt_D
        self.u_uSt = u_uSt
        self.flowType = flowType
        self._preciscion = 4

        # Pipe parameters
        self.chokedLength = nan
        self.pipeLength = nan
        self.pipeDiameter = nan
        self.frictionCoeff = nan

        # Down stream conditions if pipe parameters are given
        self.dwnStrmMach = nan
        self.dwnStrm_t_tSt = nan
        self.dwnStrm_p_pSt = nan
        self.dwnStrm_rho_rhoSt = nan
        self.dwnStrm_po_poSt = nan
        self.dwnStrm_f4LSt_D = nan
        self.dwnStrm_u_uSt = nan

        # Downstream / Initial Conditions
        self.t2_t1 = nan
        self.p2_p1 = nan
        self.rho2_rho1 = nan
        self.po2_po1 = nan
        self.f4LD2_f4LD1 = nan
        self.u2_u1 = nan

        if self.__checkValue(self.t_tSt):
            self.mach = FannoFlowRelations.calcMachFrom_T_TSt(self.t_tSt, self.gamma)

        elif self.__checkValue(self.p_pSt):
            self.mach = FannoFlowRelations.calcMachFrom_P_PSt(self.p_pSt, self.gamma)

        elif self.__checkValue(self.rho_rhoSt):
            self.mach = FannoFlowRelations.calcMachFrom_Rho_RhoSt(self.rho_rhoSt, self.gamma)

        elif self.__checkValue(self.po_poSt):
            self.mach = FannoFlowRelations.calcMachFrom_Po_PoSt(
                self.po_poSt, self.gamma, flowType=self.flowType
            )

        elif self.__checkValue(self.f4LSt_D):
            self.mach = FannoFlowRelations.calcMachFrom_4FLSt_D(
                self.f4LSt_D, self.gamma, flowType=self.flowType
            )

        elif self.__checkValue(self.u_uSt):
            self.mach = FannoFlowRelations.calcMachFrom_U_USt(self.u_uSt, self.gamma)

        if self.__checkValue(self.mach):
            self.__calculateState()

    def __str__(self) -> str:
        gammaStr = str(round(self.gamma, self._preciscion))
        machStr = str(round(self.mach, self._preciscion))
        tempStr = str(round(self.t_tSt, self._preciscion))
        prStr = str(round(self.p_pSt, self._preciscion))
        porStr = str(round(self.po_poSt, self._preciscion))
        rhoStr = str(round(self.rho_rhoSt, self._preciscion))
        f4ldStr = str(round(self.f4LSt_D, self._preciscion))
        UStr = str(round(self.u_uSt, self._preciscion))

        # Pipe Condtion Str
        chkLenStr = str(round(self.chokedLength, self._preciscion))
        pipeLenStr = str(round(self.pipeLength, self._preciscion))
        diamStr = str(round(self.pipeDiameter, self._preciscion))
        coeffStr = str(round(self.frictionCoeff, self._preciscion))

        # Final Condition Strings
        finalmachStr = str(round(self.dwnStrmMach, self._preciscion))
        finaltempStr = str(round(self.dwnStrm_t_tSt, self._preciscion))
        finalprStr = str(round(self.dwnStrm_p_pSt, self._preciscion))
        finalporStr = str(round(self.dwnStrm_po_poSt, self._preciscion))
        finalrhoStr = str(round(self.dwnStrm_rho_rhoSt, self._preciscion))
        finalf4ldStr = str(round(self.dwnStrm_f4LSt_D, self._preciscion))
        finalUStr = str(round(self.dwnStrm_u_uSt, self._preciscion))

        # Friction Jump Ratios
        jumptempStr = str(round(self.t2_t1, self._preciscion))
        jumpprStr = str(round(self.p2_p1, self._preciscion))
        jumpporStr = str(round(self.po2_po1, self._preciscion))
        jumprhoStr = str(round(self.rho2_rho1, self._preciscion))
        jumpf4ldStr = str(round(self.f4LD2_f4LD1, self._preciscion))
        jumpUStr = str(round(self.u2_u1, self._preciscion))

        width = 50 - 2  # Width minus edges
        sep = "|{:{width}}|\n".format("", width=width)

        # Initial Conditions
        header = "|{:=^{width}}|\n".format("Fanno Relations at Mach: {}".format(machStr), width=width)
        gamma = "|{:<{width}}{}|\n".format("\u03B3", gammaStr, width=width - len(gammaStr))
        tempR = "|{:-<{width}}{}|\n".format("T/T*", tempStr, width=width - len(tempStr))
        pressR = "|{:<{width}}{}|\n".format("P/P*", prStr, width=width - len(prStr))
        poR = "|{:-<{width}}{}|\n".format("Po/Po*", porStr, width=width - len(porStr))
        rhoR = "|{:<{width}}{}|\n".format("\u03C1/\u03C1*", rhoStr, width=width - len(rhoStr))
        lenR = "|{:-<{width}}{}|\n".format("4fL*/D", f4ldStr, width=width - len(f4ldStr))
        uStR = "|{:<{width}}{}|\n".format("U/U*", UStr, width=width - len(UStr))
        flow = "|{:-<{width}}{}|\n".format("FlowType", self.flowType, width=width - len(self.flowType))

        # Pipe Parameters
        paramSep = "|{:-^{width}}|\n".format("Pipe Parameters", width=width)
        if self.pipeLength > self.chokedLength:
            chokedLen = (
                Back.YELLOW
                + Fore.BLACK
                + "|{:<{width}}{}|".format("Length For Chocked Flow", chkLenStr, width=width - len(chkLenStr))
                + Style.RESET_ALL
                + "\n"
            )
        else:
            chokedLen = "|{:<{width}}{}|\n".format(
                "Length For Chocked Flow", chkLenStr, width=width - len(chkLenStr)
            )

        pipeLen = "|{:-<{width}}{}|\n".format("Pipe Length", pipeLenStr, width=width - len(pipeLenStr))
        pipeDiam = "|{:<{width}}{}|\n".format("Pipe Diameter", diamStr, width=width - len(diamStr))
        fricCoeff = "|{:-<{width}}{}|\n".format("Friction Coefficient", coeffStr, width=width - len(coeffStr))

        # Down Stream Conditions
        dwnStrmSep = "|{:-^{width}}|\n".format("Down Stream Conditions", width=width)
        dwnmachStr = "|{:<{width}}{}|\n".format("Mach", finalmachStr, width=width - len(finalmachStr))
        dwntempR = "|{:-<{width}}{}|\n".format("T/T*", finaltempStr, width=width - len(finaltempStr))
        dwnpressR = "|{:<{width}}{}|\n".format("P/P*", finalprStr, width=width - len(finalprStr))
        dwnpoR = "|{:-<{width}}{}|\n".format("Po/Po*", finalporStr, width=width - len(finalporStr))
        dwnrhoR = "|{:<{width}}{}|\n".format("\u03C1/\u03C1*", finalrhoStr, width=width - len(finalrhoStr))
        dwnlenR = "|{:-<{width}}{}|\n".format("4fL*/D", finalf4ldStr, width=width - len(finalf4ldStr))
        dwnuStR = "|{:<{width}}{}|\n".format("U/U*", finalUStr, width=width - len(finalUStr))

        # Parameter Ratios
        jumpSep = "|{:-^{width}}|\n".format("Conditions Across Friciton Area", width=width)
        jumppRatio = "|{:<{width}}{}|\n".format("p2/p1:", jumpprStr, width=width - len(jumpprStr))
        jumprhoRatio = "|{:-<{width}}{}|\n".format(
            "\u03C12/\u03C11:", jumprhoStr, width=width - len(jumprhoStr)
        )
        jumptempRatio = "|{:<{width}}{}|\n".format("T2/T1:", jumptempStr, width=width - len(jumptempStr))
        jumppoRatio = "|{:-<{width}}{}|\n".format("po2/po1:", jumpporStr, width=width - len(jumpporStr))
        jumpFrRatio = "|{:<{width}}{}|\n".format(
            "(4fL*/D)2 / (4fL*/D)1:", jumpf4ldStr, width=width - len(jumpf4ldStr)
        )
        jumpVelRatio = "|{:<{width}}{}|\n".format("U2/U1:", jumpUStr, width=width - len(jumpUStr))

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
                lenR,
                uStR,
                flow,
                sep,
                paramSep,
                chokedLen,
                pipeLen,
                pipeDiam,
                fricCoeff,
                sep,
                dwnStrmSep,
                dwnmachStr,
                dwntempR,
                dwnpressR,
                dwnpoR,
                dwnrhoR,
                dwnlenR,
                dwnuStR,
                sep,
                jumpSep,
                jumppRatio,
                jumprhoRatio,
                jumptempRatio,
                jumppoRatio,
                jumpFrRatio,
                jumpVelRatio,
                sep,
                finish,
            ]
        )

    def applyPipeParameters(self, diameter: float, length: float, frictionCoeff: float = 0.005) -> None:
        self.pipeDiameter = diameter
        self.pipeLength = length
        self.frictionCoeff = frictionCoeff
        self.__calculateDownStreamState()

    def __calculateState(self) -> None:
        self.t_tSt = FannoFlowRelations.calcT_Tstar(self.mach, self.gamma)
        self.p_pSt = FannoFlowRelations.calcP_Pstar(self.mach, self.gamma)
        self.rho_rhoSt = FannoFlowRelations.calcRho_RhoStar(self.mach, self.gamma)
        self.po_poSt = FannoFlowRelations.calcRho_RhoStar(self.mach, self.gamma)
        self.f4LSt_D = FannoFlowRelations.calc4FLSt_D(self.mach, self.gamma)
        self.u_uSt = FannoFlowRelations.calcU_Ustar(self.mach, self.gamma)

        if self.mach < 1:
            self.flowType = "Subsonic"
        else:
            self.flowType = "Supersonic"

    def __calculateDownStreamState(self) -> None:
        if (
            not self.__checkValue(self.pipeDiameter)
            or not self.__checkValue(self.pipeLength)
            or not self.__checkValue(self.frictionCoeff)
        ):
            return

        self.chokedLength = self.f4LSt_D * self.pipeDiameter / 4 / self.frictionCoeff
        self.dwnStrm_f4LSt_D = self.f4LSt_D - 4 * self.frictionCoeff * self.pipeLength / self.pipeDiameter

        if self.pipeLength > self.chokedLength:
            self.dwnStrmMach = 1
        else:
            self.dwnStrmMach = FannoFlowRelations.calcMachFrom_4FLSt_D(
                self.dwnStrm_f4LSt_D, self.gamma, self.flowType
            )

        self.dwnStrm_t_tSt = FannoFlowRelations.calcT_Tstar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_p_pSt = FannoFlowRelations.calcP_Pstar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_rho_rhoSt = FannoFlowRelations.calcRho_RhoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_po_poSt = FannoFlowRelations.calcRho_RhoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_u_uSt = FannoFlowRelations.calcU_Ustar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_f4LSt_D = FannoFlowRelations.calc4FLSt_D(self.dwnStrmMach, self.gamma)

        # Calculate parameter ratios across the condition
        self.t2_t1 = self.dwnStrm_t_tSt / self.t_tSt
        self.p2_p1 = self.dwnStrm_p_pSt / self.p_pSt
        self.rho2_rho1 = self.dwnStrm_rho_rhoSt / self.rho_rhoSt
        self.po2_po1 = self.dwnStrm_po_poSt / self.po_poSt
        self.f4LD2_f4LD1 = self.dwnStrm_f4LSt_D / self.f4LSt_D
        self.u2_u1 = self.dwnStrm_u_uSt / self.u_uSt

    def __checkValue(self, var: float) -> bool:
        if isnan(var):
            return False

        if var < 0:
            return False

        return True

    @staticmethod
    def calcT_Tstar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates T_T* given gamma and Mach, offset can be applied for root finding"""
        return (gamma + 1) / (2 + (gamma - 1) * pow(mach, 2)) - offset

    @staticmethod
    def calcMachFrom_T_TSt(t_tSt: float, gamma: float) -> float:
        """Calculates Mach given a T_T* value and gamma"""
        return brenth(FannoFlowRelations.calcT_Tstar, 1e-9, 40, args=(gamma, t_tSt,))

    @staticmethod
    def calcP_Pstar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates p/p* given gamma and mach"""
        return sqrt(FannoFlowRelations.calcT_Tstar(mach, gamma)) / mach - offset

    @staticmethod
    def calcMachFrom_P_PSt(p_pSt: float, gamma: float) -> float:
        """ Calculates mach given p/p* and gamma"""
        return brenth(FannoFlowRelations.calcP_Pstar, 1e-9, 40, args=(gamma, p_pSt,))

    @staticmethod
    def calcRho_RhoStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates rho/rho* given gamma and mach"""
        return sqrt(1 / FannoFlowRelations.calcT_Tstar(mach, gamma)) / mach - offset

    @staticmethod
    def calcMachFrom_Rho_RhoSt(rho_rhoSt: float, gamma: float) -> float:
        """ Calculates mach given rho/rho* and gamma"""
        return brenth(FannoFlowRelations.calcRho_RhoStar, 1e-9, 40, args=(gamma, rho_rhoSt,),)

    @staticmethod
    def calcRho_RhoStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates po/po* given mach and gamma"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        return pow(1 / FannoFlowRelations.calcT_Tstar(mach, gamma), gp1 / (2 * gm1)) / mach - offset

    @staticmethod
    def calcMachFrom_Po_PoSt(po_poSt: float, gamma: float, flowType: str = "Supersonic") -> float:
        """ Calculates mach given po/po* and gamma"""
        tolerance = 1e-5
        if po_poSt == 1.0:
            return 1
        elif flowType == "Supersonic":
            return brenth(FannoFlowRelations.calcRho_RhoStar, 1 + tolerance, 40, args=(gamma, po_poSt,),)
        elif flowType == "Subsonic":
            return brenth(
                FannoFlowRelations.calcRho_RhoStar, tolerance, 1 - tolerance, args=(gamma, po_poSt,),
            )
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType) + Back.RESET + Fore.RESET
            )

    @staticmethod
    def calc4FLSt_D(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates 4fL*/D given mach and gamma"""
        gp1 = gamma + 1
        t_tSt = FannoFlowRelations.calcT_Tstar(mach, gamma)
        mSqr = pow(mach, 2)
        return (1 - mSqr) / (gamma * mSqr) + (gp1 / (2 * gamma)) * log(t_tSt * mSqr) - offset

    @staticmethod
    def calcMachFrom_4FLSt_D(f4lSt_d: float, gamma: float, flowType: str = "Supersonic") -> float:
        """ Calculates mach given gamma and 4FL*/D"""
        if f4lSt_d == 0.0:
            return 1
        elif flowType == "Supersonic":
            return brenth(FannoFlowRelations.calc4FLSt_D, 1.00001, 50, args=(gamma, f4lSt_d,),)
        elif flowType == "Subsonic":
            return brenth(FannoFlowRelations.calc4FLSt_D, 1e-5, 0.9999, args=(gamma, f4lSt_d,),)
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType) + Back.RESET + Fore.RESET
            )

    @staticmethod
    def calcU_Ustar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates U_U* given mach and gamma"""
        t_tSt = FannoFlowRelations.calcT_Tstar(mach, gamma)
        return mach * sqrt(t_tSt) - offset

    @staticmethod
    def calcMachFrom_U_USt(u_uSt: float, gamma: float) -> float:
        """ Calculates Mach given U_U* and gamma"""
        return brenth(FannoFlowRelations.calcU_Ustar, 1e-9, 40, args=(gamma, u_uSt))

