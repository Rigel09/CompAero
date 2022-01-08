from CompAero.internal import (
    InvalidOptionCombinationError,
    checkValue,
    FlowState,
    footer,
    named_header,
    seperator,
    to_string,
    GammaNotDefinedError,
)
from CompAero.greek_letters import LowerCaseGreek as lcg

from math import sqrt, nan, pow
from scipy.optimize import brenth


class IsentropicRelations:
    """ This class is a collective name space for basic calculations regarding isentropic flows. 
    The constructor of this class can also determine the entire state of the flow given a partial state of the flow 

    Args:
        gamma (float): ratio of specific heats
        mach (float, optional): Mach number of the flow. Defaults to nan.
        p0_p (float, optional): Ratio of total pressure to static pressure. Defaults to nan.
        t0_t (float, optional): Ratio of total temperature to static temperature. Defaults to nan.
        rho0_rho (float, optional): Ratio of total density to static density. Defaults to nan.
        a_aStar (float, optional): Ratio of Nozzle Area to sonic throat area. Defaults to nan.
        flowType (FlowState, optional): States wether the flow is subsonic or supersonic. Used for Area Ratio Calculations. Defaults to FlowState.SUPER_SONIC.

    Raises:
        InvalidOptionCombinationError: Raised if invalid combination of values are given to the constructor
        
    Useage:
        To use this class pass gamma and one of the known parameters of the flow and the rest are calculated. 
    
    Valid Combinations of Parameters:
        gamma, mach
        gamma, p0/p
        gamma, t0/t
        gamma, rho0/rho
        gamma, a/a*, flowtype (flow type defaults to super sonic)
    """

    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        p0_p: float = nan,
        t0_t: float = nan,
        rho0_rho: float = nan,
        a_aStar: float = nan,
        flowType: FlowState = FlowState.SUPER_SONIC,
    ) -> None:
        self.gamma = gamma
        self.mach = mach
        self.p0_p = p0_p
        self.t0_t = t0_t
        self.rho0_rho = rho0_rho
        self.a_aStar = a_aStar
        self.flowType = flowType
        self.precision = 4

        # Calculate parameters based on what was passed in
        if not checkValue(self.gamma):
            raise GammaNotDefinedError()

        if checkValue(self.mach):
            pass
        elif checkValue(self.p0_p):
            self.mach = IsentropicRelations.calc_mach_from_p0_p(self.p0_p, self.gamma)
        elif checkValue(self.t0_t):
            self.mach = IsentropicRelations.calc_mach_from_T0_T(self.t0_t, self.gamma)
        elif checkValue(self.rho0_rho):
            self.mach = IsentropicRelations.calc_mach_from_rho0_rho(self.rho0_rho, self.gamma)
        elif checkValue(self.a_aStar):
            self.mach = IsentropicRelations.calc_mach_from_A_Astar(self.a_aStar, self.gamma, self.flowType)
        else:
            raise InvalidOptionCombinationError()

        if checkValue(self.mach):
            self.__calculateStateFromMach()

    def __calculateStateFromMach(self) -> None:
        self.t0_t = IsentropicRelations.calc_T0_T(self.mach, self.gamma)
        self.p0_p = IsentropicRelations.calc_P0_P(self.mach, self.gamma)
        self.rho0_rho = IsentropicRelations.calc_rho0_rho(self.mach, self.gamma)
        self.a_aStar = IsentropicRelations.calc_A_Astar(self.mach, self.gamma)
        self.flowType = FlowState.SUPER_SONIC if self.mach > 1.0 else FlowState.SUB_SONIC

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Isentropic Flow State at Mach", self.mach, self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision, dot_line=True),
                to_string("p0/p", self.p0_p, self.precision),
                to_string("T0/T", self.t0_t, self.precision, dot_line=True),
                to_string("{}0/{}".format(lcg.rho, lcg.rho), self.rho0_rho, self.precision),
                to_string("A/A*", self.a_aStar, self.precision, dot_line=True),
                footer(),
            ]
        )

    @staticmethod
    def calc_T0_T(mach: float, gamma: float, ofFlowStateet: float = 0.0) -> float:
        """Calculates the ratio of Total Temperature to static temperature (T0/T)

        Args:
            mach (float): Mach number of the flow
            gamma (float): Ratio of specific heats of the flow
            ofFlowStateet (float, optional): OfFlowStateet that can be used for root finding for a specific value. Defaults to 0.0.

        Returns:
            float: T0/T
        """

        return 1 + (gamma - 1) / 2 * pow(mach, 2) - ofFlowStateet

    @staticmethod
    def calc_mach_from_T0_T(t0_t: float, gamma: float) -> float:
        """Calculates the Mach number for a flow given the ratio of total temperature to static temperature

        Args:
            t0_t (float): Ratio of total temperature to static temperature
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return sqrt((t0_t - 1) * 2 / (gamma - 1))

    @staticmethod
    def calc_P0_P(mach: float, gamma: float, ofFlowStateet: float = 0.0) -> float:
        """Calculates the ratio of Total pressure to static pressure (P0/P)

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            ofFlowStateet (float, optional): [description]. Defaults to 0.0.

        Returns:
            float: P0/P
        """

        return pow((1 + (gamma - 1) / 2 * pow(mach, 2)), gamma / (gamma - 1)) - ofFlowStateet

    @staticmethod
    def calc_mach_from_p0_p(p0_p: float, gamma: float) -> float:
        """Calculates the Mach number for a flow given the ratio of total pressure to static pressure

        Args:
            p0_p (float): Ratio of total pressure to static pressure
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """

        return brenth(IsentropicRelations.calc_P0_P, 0, 30, args=(gamma, p0_p))

    @staticmethod
    def calc_rho0_rho(mach: float, gamma: float, ofFlowStateet: float = 0.0) -> float:
        """Calculates the ratio of Total density to static density (rho0/rho)

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            ofFlowStateet (float, optional): [description]. Defaults to 0.0.

        Returns:
            float: rho0/rho
        """

        return pow((1 + (gamma - 1) / 2 * pow(mach, 2)), 1 / (gamma - 1)) - ofFlowStateet

    @staticmethod
    def calc_mach_from_rho0_rho(rho0_rho: float, gamma: float) -> float:
        """Calculates the Mach number for a flow given the ratio of total density to static density

        Args:
            rho0_rho (float): Ratio of total density to static density
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """

        return brenth(IsentropicRelations.calc_rho0_rho, 0, 30, args=(gamma, rho0_rho))

    @staticmethod
    def calc_A_Astar(mach: float, gamma: float, ofFlowStateet: float = 0.0) -> float:
        """Calculates the ratio of Nozzle Area to Sonic Throat Area (A/A*)

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            ofFlowStateet (float, optional): [description]. Defaults to 0.0.

        Returns:
            float: A/A*
        """

        gm1 = gamma - 1
        gp1 = gamma + 1
        mSqr = pow(mach, 2)

        nonRaised = 2 / gp1 * (1 + gm1 / 2 * mSqr)

        return sqrt(pow(nonRaised, gp1 / gm1) / mSqr) - ofFlowStateet

    @staticmethod
    def calc_mach_from_A_Astar(
        A_Astar: float, gamma: float, flowType: FlowState = FlowState.SUPER_SONIC
    ) -> float:
        """Calculates the mach number for a flow given the nozzle area ratio and flow type

        Args:
            A_Astar (float): Ratio of nozzle area ratio to sonic area ratio
            gamma (float): ratio of specific heats
            flowType (FlowState, optional): Type of flow whether it is super sonic of subsonic. Defaults to FlowState.SUPER_SONIC.

        Returns:
            float: mach number
        """
        assert isinstance(flowType, FlowState)
        if A_Astar == 1.0:
            return A_Astar
        elif flowType == FlowState.SUPER_SONIC:
            return brenth(IsentropicRelations.calc_A_Astar, 1, 30, args=(gamma, A_Astar))
        elif flowType == FlowState.SUB_SONIC:
            return brenth(IsentropicRelations.calc_A_Astar, 0.001, 1, args=(gamma, A_Astar))
        else:
            print("Unsupported flow type passed to A/A* calculations. Type: {}".format(flowType.name))
            return nan
