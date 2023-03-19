from math import sqrt, nan, pow, isnan, log
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from enum import Enum
from CompAero.greek_letters import LowerCaseGreek as lcg
from CompAero.internal import (
    FlowState,
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    checkValue,
    to_string,
    footer,
    named_header,
    named_subheader,
    seperator,
)

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
    """ This class is a collective name space for basic calculations regarding Fanno flows. 
        The constructor of this class can also determine the entire state of the flow given a partial state of the flow 

    Args:
        gamma (float): ratio of specific heats      
        mach (float, optional): mach number of the flow. Defaults to nan.
        t_tSt (float, optional): Ratio of Temperature to sonic temperature. Defaults to nan.
        p_pSt (float, optional): Ratio of Pressure to sonic pressure. Defaults to nan.
        rho_rhoSt (float, optional): Ratio of density to sonic density. Defaults to nan.
        po_poSt (float, optional): Ratio of Total pressure to sonic total pressre. Defaults to nan.
        f4LSt_D (float, optional): Effect of friction on pipe. Defaults to nan.
        u_uSt (float, optional): Velocity to sonic Velocity. Defaults to nan.
        flowType (FlowState, optional):  States wether the flow is subsonic or supersonic. Used for Area Ratio Calculations. Defaults to FlowState.SUPER_SONIC.
        
    Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and flow state cannot be determined
        
    Useage:
        To use this class pass gamma and one of the known parameters of the flow and the rest are calculated. 

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
        t_tSt: float = nan,
        p_pSt: float = nan,
        rho_rhoSt: float = nan,
        po_poSt: float = nan,
        f4LSt_D: float = nan,
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
        self.f4LSt_D = f4LSt_D
        """ Friction parameter 4fL*/D"""
        self.u_uSt = u_uSt
        """ Ratio of velocity to the sonic velcoity U/U*"""
        self.flowType = flowType
        """ Type of flow which is either subsonic or supersonic (Type: flowstate) """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        # Pipe parameters
        self.chokedLength = nan
        """ Pipe lenght required to choke the flow"""
        self.pipeLength = nan
        """ Current pipe lenght"""
        self.pipeDiameter = nan
        """ Diameter of the pipe"""
        self.frictionCoeff = nan
        """ Friction coefficient"""

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
        self.dwnStrm_f4LSt_D = nan
        """ Friction parameter 4fL*/D of the downstream flow"""
        self.dwnStrm_u_uSt = nan
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
        self.f4LD2_f4LD1 = nan
        """ Ratio of the downstream friction parameter to the upstream friction parameter"""
        self.u2_u1 = nan
        """ Ratio of downstream velocity to upstream velocity """

        if not checkValue(self.gamma):
            raise GammaNotDefinedError()

        if checkValue(self.mach):
            pass

        elif checkValue(self.t_tSt):
            self.mach = FannoFlowRelations.calc_mach_from_T_TStar(self.t_tSt, self.gamma)

        elif checkValue(self.p_pSt):
            self.mach = FannoFlowRelations.calc_mach_from_P_PStar(self.p_pSt, self.gamma)

        elif checkValue(self.rho_rhoSt):
            self.mach = FannoFlowRelations.calc_mach_from_Rho_RhoStar(self.rho_rhoSt, self.gamma)

        elif checkValue(self.po_poSt):
            self.mach = FannoFlowRelations.calc_mach_from_Po_PoStar(
                self.po_poSt, self.gamma, flowType=self.flowType
            )

        elif checkValue([self.f4LSt_D, self.flowType]):
            self.mach = FannoFlowRelations.calc_mach_from_4FLSt_D(
                self.f4LSt_D, self.gamma, flowType=self.flowType
            )

        elif checkValue(self.u_uSt):
            self.mach = FannoFlowRelations.calc_mach_from_U_USt(self.u_uSt, self.gamma)

        else:
            raise InvalidOptionCombinationError()

        if checkValue(self.mach):
            self.__calculateState()

    @property
    def chockedFlow(self) -> bool:
        """ True if the flow has reached the choked condition (I.E. Sonic) """
        return self.pipeLength >= self.chokedLength

    def __str__(self) -> str:
        color = Back.GREEN + Fore.BLACK if self.pipeLength < self.chokedLength else Back.YELLOW + Fore.BLACK

        return "".join(
            [
                named_header("Fanno Relations at Mach", self.mach, precision=self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string("T/T*", self.t_tSt, self.precision, dot_line=True),
                to_string("P/P*", self.p_pSt, self.precision),
                to_string("{}/{}*".format(lcg.rho, lcg.rho), self.rho_rhoSt, self.precision, dot_line=True),
                to_string("4FL*/D", self.f4LSt_D, self.precision),
                to_string("U/U*", self.u_uSt, self.precision, dot_line=True),
                to_string("Flow Type", self.flowType.name, self.precision),
                seperator(),
                named_subheader("Pipe Parameters"),
                to_string("Length For Chocked Flow", self.chokedLength, self.precision),
                color,
                to_string("Is Flow Choked? ", self.chockedFlow, self.precision, dot_line=True),
                to_string("Pipe Length", self.pipeLength, self.precision),
                to_string("Pipe Diameter", self.pipeDiameter, self.precision, dot_line=True),
                to_string("Friction Coefficient", self.frictionCoeff, self.precision),
                seperator(),
                named_subheader("Down Stream Conditions"),
                to_string("Mach", self.dwnStrmMach, self.precision),
                to_string("T/T*", self.t_tSt, self.precision, dot_line=True),
                to_string("P/P*", self.dwnStrm_p_pSt, self.precision),
                to_string("P0/P0*", self.dwnStrm_po_poSt, self.precision, dot_line=True),
                to_string("{}/{}*".format(lcg.rho, lcg.rho), self.dwnStrm_rho_rhoSt, self.precision),
                to_string("4FL*/D", self.dwnStrm_f4LSt_D, self.precision, dot_line=True),
                to_string("U/U*", self.dwnStrm_u_uSt, self.precision),
                seperator(),
                named_subheader("Conditions Across Friction Area"),
                to_string("p2/p1", self.p2_p1, self.precision),
                to_string("{}2/{}1".format(lcg.rho, lcg.rho), self.rho2_rho1, self.precision, dot_line=True),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("4FL*/D2 / 4FL*/D 1", self.f4LD2_f4LD1, self.precision),
                to_string("U2/U1", self.u2_u1, self.precision, dot_line=True),
                footer(),
            ]
        )

    def apply_pipe_parameters(self, diameter: float, length: float, frictionCoeff: float = 0.005) -> None:
        """This functions applies parameters of a known pipe to the determined state of the flow. 
            This allows the state at the downstream end of the pipe or pipe section to be found

        Args:
            diameter (float): Diameter of pipe
            length (float): Length of pipe
            frictionCoeff (float, optional): Friction coefficient of the pipe. Defaults to 0.005 which holds for Reynolds > 10^5.
        """
        self.pipeDiameter = diameter
        self.pipeLength = length
        self.frictionCoeff = frictionCoeff
        self.__calculateDownStreamState()

    def __calculateState(self) -> None:
        self.t_tSt = FannoFlowRelations.calc_T_Tstar(self.mach, self.gamma)
        self.p_pSt = FannoFlowRelations.calc_P_Pstar(self.mach, self.gamma)
        self.rho_rhoSt = FannoFlowRelations.calc_Rho_RhoStar(self.mach, self.gamma)
        self.po_poSt = FannoFlowRelations.calc_Po_PoStar(self.mach, self.gamma)
        self.f4LSt_D = FannoFlowRelations.calc_4FLSt_D(self.mach, self.gamma)
        self.u_uSt = FannoFlowRelations.calc_U_UStar(self.mach, self.gamma)

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
            self.dwnStrmMach = FannoFlowRelations.calc_mach_from_4FLSt_D(
                self.dwnStrm_f4LSt_D, self.gamma, self.flowType
            )

        self.dwnStrm_t_tSt = FannoFlowRelations.calc_T_Tstar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_p_pSt = FannoFlowRelations.calc_P_Pstar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_rho_rhoSt = FannoFlowRelations.calc_Rho_RhoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_po_poSt = FannoFlowRelations.calc_Po_PoStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_u_uSt = FannoFlowRelations.calc_U_UStar(self.dwnStrmMach, self.gamma)
        self.dwnStrm_f4LSt_D = FannoFlowRelations.calc_4FLSt_D(self.dwnStrmMach, self.gamma)

        # Calculate parameter ratios across the condition
        self.t2_t1 = self.dwnStrm_t_tSt / self.t_tSt
        self.p2_p1 = self.dwnStrm_p_pSt / self.p_pSt
        self.rho2_rho1 = self.dwnStrm_rho_rhoSt / self.rho_rhoSt
        self.po2_po1 = self.dwnStrm_po_poSt / self.po_poSt
        self.f4LD2_f4LD1 = self.dwnStrm_f4LSt_D / self.f4LSt_D
        self.u2_u1 = self.dwnStrm_u_uSt / self.u_uSt

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
        return (gamma + 1) / (2 + (gamma - 1) * pow(mach, 2)) - offset

    @staticmethod
    def calc_mach_from_T_TStar(t_tSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of static temperature to sonic static temperature T/T*

        Args:
            t_tSt (float): Ratio of static temperature to sonic static temperature T/T*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(FannoFlowRelations.calc_T_Tstar, 1e-9, 40, args=(gamma, t_tSt,))

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
        return sqrt(FannoFlowRelations.calc_T_Tstar(mach, gamma)) / mach - offset

    @staticmethod
    def calc_mach_from_P_PStar(p_pSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of static pressure to sonic static pressure P/P*

        Args:
            p_pSt (float): Ratio of static pressure to sonic static pressure P/P*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(FannoFlowRelations.calc_P_Pstar, 1e-9, 40, args=(gamma, p_pSt,))

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
        return sqrt(1 / FannoFlowRelations.calc_T_Tstar(mach, gamma)) / mach - offset

    @staticmethod
    def calc_mach_from_Rho_RhoStar(rho_rhoSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of density to sonic density Rho/Rho*

        Args:
            rho_rhoSt (float): Ratio of density to sonic density Rho/Rho*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(FannoFlowRelations.calc_Rho_RhoStar, 1e-9, 40, args=(gamma, rho_rhoSt,),)

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
        return pow(1 / FannoFlowRelations.calc_T_Tstar(mach, gamma), gp1 / (2 * gm1)) / mach - offset

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
            return brenth(FannoFlowRelations.calc_Po_PoStar, 1 + tolerance, 40, args=(gamma, po_poSt,),)
        elif flowType == FlowState.SUB_SONIC:
            return brenth(
                FannoFlowRelations.calc_Po_PoStar, tolerance, 1 - tolerance, args=(gamma, po_poSt,),
            )
        else:
            raise ValueError(
                Back.RED + Fore.BLACK + "Flow Type [{}] not supported for Fanno"
                " Flow! Accepted Types: 'Supersonic' & 'Subsonic'".format(flowType.name)
                + Back.RESET
                + Fore.RESET
            )

    @staticmethod
    def calc_4FLSt_D(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates friction parameter for flow

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: 4FL*/D
        """
        gp1 = gamma + 1
        t_tSt = FannoFlowRelations.calc_T_Tstar(mach, gamma)
        mSqr = pow(mach, 2)
        return (1 - mSqr) / (gamma * mSqr) + (gp1 / (2 * gamma)) * log(t_tSt * mSqr) - offset

    @staticmethod
    def calc_mach_from_4FLSt_D(
        f4lSt_d: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """ Calculates the mach number from the friction parameter

        Args:
            f4lSt_d (float): friction parameter 4FL*/D
            gamma (float): ratio of specific heats
            flowType (FlowState, optional):  Type of flow whether it is super sonic of subsonic. Defaults to FlowState.SUPER_SONIC.

        Raises:
            ValueError: Raised if Flow State is not supersonic of subsonic  

        Returns:
            float: mach number
        """
        if f4lSt_d == 0.0:
            return 1
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(FannoFlowRelations.calc_4FLSt_D, 1.00001, 50, args=(gamma, f4lSt_d,),)
        elif flowType == FlowState.SUB_SONIC:
            return brenth(FannoFlowRelations.calc_4FLSt_D, 1e-5, 0.9999, args=(gamma, f4lSt_d,),)
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
        t_tSt = FannoFlowRelations.calc_T_Tstar(mach, gamma)
        return mach * sqrt(t_tSt) - offset

    @staticmethod
    def calc_mach_from_U_USt(u_uSt: float, gamma: float) -> float:
        """Calculates the mach number based of the ratio of velocity to sonic velocity U/U*

        Args:
            u_uSt (float): Ratio of velocity to sonic velocity U/U*
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(FannoFlowRelations.calc_U_UStar, 1e-9, 40, args=(gamma, u_uSt))

