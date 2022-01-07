from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from CompAero.greek_letters import LowerCaseGreek as lcg
from CompAero.internal import (
    FlowState,
    checkValue,
    data_value_to_string,
    footer,
    named_header,
    named_subheader,
    seperator,
)


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
        flowType: FlowState = FlowState.SUPER_SONIC,
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

        if checkValue(self.t_tSt):
            self.mach = FannoFlowRelations.calcMachFrom_T_TSt(self.t_tSt, self.gamma)

        elif checkValue(self.p_pSt):
            self.mach = FannoFlowRelations.calcMachFrom_P_PSt(self.p_pSt, self.gamma)

        elif checkValue(self.rho_rhoSt):
            self.mach = FannoFlowRelations.calcMachFrom_Rho_RhoSt(self.rho_rhoSt, self.gamma)

        elif checkValue(self.po_poSt):
            self.mach = FannoFlowRelations.calcMachFrom_Po_PoSt(
                self.po_poSt, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.f4LSt_D):
            self.mach = FannoFlowRelations.calcMachFrom_4FLSt_D(
                self.f4LSt_D, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.u_uSt):
            self.mach = FannoFlowRelations.calcMachFrom_U_USt(self.u_uSt, self.gamma)

        if checkValue(self.mach):
            self.__calculateState()

    @property
    def chockedFlow(self) -> bool:
        return self.pipeLength >= self.chokedLength

    def __str__(self) -> str:
        color = "" if self.pipeLength < self.chokedLength else Back.YELLOW + Fore.BLACK

        return "".join(
            [
                named_header("Fanno Relations at Mach", self.mach, precision=self._preciscion),
                seperator(),
                data_value_to_string("{}:".format(lcg.gamma), self.gamma, self._preciscion),
                data_value_to_string("T/T*:", self.t_tSt, self._preciscion, dot_line=True),
                data_value_to_string("P/P*:", self.p_pSt, self._preciscion),
                data_value_to_string(
                    "{}/{}*:".format(lcg.rho, lcg.rho), self.rho_rhoSt, self._preciscion, dot_line=True
                ),
                data_value_to_string("4FL*/D:", self.f4LSt_D, self._preciscion),
                data_value_to_string("U/U*:", self.u_uSt, self._preciscion, dot_line=True),
                data_value_to_string("Flow Type:", self.flowType.name, self._preciscion),
                seperator(),
                named_subheader("Pipe Parameters"),
                color,
                data_value_to_string("Length For Chocked Flow:", self.chokedLength, self._preciscion),
                data_value_to_string("Is Flow Choked? ", self.chockedFlow, self._preciscion, dot_line=True),
                data_value_to_string("Pipe Length:", self.pipeLength, self._preciscion),
                data_value_to_string("Pipe Diameter:", self.pipeDiameter, self._preciscion, dot_line=True),
                data_value_to_string("Friction Coefficient:", self.frictionCoeff, self._preciscion),
                seperator(),
                named_subheader("Down Stream Conditions"),
                data_value_to_string("Mach:", self.dwnStrmMach, self._preciscion),
                data_value_to_string("T/T*:", self.t_tSt, self._preciscion, dot_line=True),
                data_value_to_string("P/P*:", self.dwnStrm_p_pSt, self._preciscion),
                data_value_to_string("P0/P0*:", self.dwnStrm_po_poSt, self._preciscion, dot_line=True),
                data_value_to_string(
                    "{}/{}*:".format(lcg.rho, lcg.rho), self.dwnStrm_rho_rhoSt, self._preciscion
                ),
                data_value_to_string("4FL*/D:", self.dwnStrm_f4LSt_D, self._preciscion, dot_line=True),
                data_value_to_string("U/U*:", self.dwnStrm_u_uSt, self._preciscion),
                seperator(),
                named_subheader("Conditions Across Friction Area"),
                data_value_to_string("p2/p1:", self.p2_p1, self._preciscion),
                data_value_to_string(
                    "{}2/{}1".format(lcg.rho, lcg.rho), self.rho2_rho1, self._preciscion, dot_line=True
                ),
                data_value_to_string("T2/T1:", self.t2_t1, self._preciscion),
                data_value_to_string("P02/P01:", self.po2_po1, self._preciscion, dot_line=True),
                data_value_to_string("4FL*/D2 / 4FL*/D 1:", self.f4LD2_f4LD1, self._preciscion),
                data_value_to_string("U2/U1:", self.u2_u1, self._preciscion, dot_line=True),
                seperator(),
                footer(),
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
        self.po_poSt = FannoFlowRelations.calcPo_PoStar(self.mach, self.gamma)
        self.f4LSt_D = FannoFlowRelations.calc4FLSt_D(self.mach, self.gamma)
        self.u_uSt = FannoFlowRelations.calcU_Ustar(self.mach, self.gamma)

        if self.mach < 1:
            self.flowType = FlowState.SUB_SONIC
        else:
            self.flowType = FlowState.SUPER_SONIC

    def __calculateDownStreamState(self) -> None:
        if not checkValue([self.pipeDiameter, self.pipeLength, self.frictionCoeff]):
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
        self.dwnStrm_po_poSt = FannoFlowRelations.calcPo_PoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_u_uSt = FannoFlowRelations.calcU_Ustar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_f4LSt_D = FannoFlowRelations.calc4FLSt_D(self.dwnStrmMach, self.gamma)

        # Calculate parameter ratios across the condition
        self.t2_t1 = self.dwnStrm_t_tSt / self.t_tSt
        self.p2_p1 = self.dwnStrm_p_pSt / self.p_pSt
        self.rho2_rho1 = self.dwnStrm_rho_rhoSt / self.rho_rhoSt
        self.po2_po1 = self.dwnStrm_po_poSt / self.po_poSt
        self.f4LD2_f4LD1 = self.dwnStrm_f4LSt_D / self.f4LSt_D
        self.u2_u1 = self.dwnStrm_u_uSt / self.u_uSt

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
    def calcPo_PoStar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates po/po* given mach and gamma"""
        gp1 = gamma + 1
        gm1 = gamma - 1
        return pow(1 / FannoFlowRelations.calcT_Tstar(mach, gamma), gp1 / (2 * gm1)) / mach - offset

    @staticmethod
    def calcMachFrom_Po_PoSt(
        po_poSt: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """ Calculates mach given po/po* and gamma"""
        tolerance = 1e-5
        if po_poSt == 1.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(FannoFlowRelations.calcPo_PoStar, 1 + tolerance, 40, args=(gamma, po_poSt,),)
        elif flowType == FlowState.SUB_SONIC:
            return brenth(FannoFlowRelations.calcPo_PoStar, tolerance, 1 - tolerance, args=(gamma, po_poSt,),)
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calc4FLSt_D(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Calculates 4fL*/D given mach and gamma"""
        gp1 = gamma + 1
        t_tSt = FannoFlowRelations.calcT_Tstar(mach, gamma)
        mSqr = pow(mach, 2)
        return (1 - mSqr) / (gamma * mSqr) + (gp1 / (2 * gamma)) * log(t_tSt * mSqr) - offset

    @staticmethod
    def calcMachFrom_4FLSt_D(
        f4lSt_d: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """ Calculates mach given gamma and 4FL*/D"""
        if f4lSt_d == 0.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(FannoFlowRelations.calc4FLSt_D, 1.00001, 50, args=(gamma, f4lSt_d,),)
        elif flowType == FlowState.SUB_SONIC:
            return brenth(FannoFlowRelations.calc4FLSt_D, 1e-5, 0.9999, args=(gamma, f4lSt_d,),)
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
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

