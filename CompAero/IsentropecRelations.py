from math import sqrt, nan, pow, isnan
import numpy as np
from scipy.optimize import brenth
from CompAero.internal import (
    checkValue,
    FlowState as FS,
    footer,
    named_header,
    named_subheader,
    seperator,
    data_value_to_string,
)
from CompAero.greek_letters import LowerCaseGreek as lcg


class IsentropicRelations:
    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        p0_p: float = nan,
        t0_t: float = nan,
        rho0_rho: float = nan,
        a_aStar: float = nan,
        flowType: FS = FS.SUPER_SONIC,
    ) -> None:
        self.gamma = gamma
        self.mach = mach
        self.p0_p = p0_p
        self.t0_t = t0_t
        self.rho0_rho = rho0_rho
        self.a_aStar = a_aStar
        self.flowType = flowType
        self._precision = 4

        # Calculate parameters based on what was passed in
        if not checkValue(self.gamma):
            return

        if checkValue(self.p0_p):
            self.mach = IsentropicRelations.calcMachFrom_p0_p(self.p0_p, self.gamma)
        elif checkValue(self.t0_t):
            self.mach = IsentropicRelations.calcMachFrom_T0_T(self.t0_t, self.gamma)
        elif checkValue(self.rho0_rho):
            self.mach = IsentropicRelations.calcMachFrom_rho0_rho(self.rho0_rho, self.gamma)
        elif checkValue(self.a_aStar):
            self.mach = IsentropicRelations.calcMachFrom_A_Astar(self.a_aStar, self.gamma, self.flowType)

        if checkValue(self.mach):
            self.__calculateStateFromMach()

    def __calculateStateFromMach(self) -> None:
        self.t0_t = IsentropicRelations.calc_T0_T(self.mach, self.gamma)
        self.p0_p = IsentropicRelations.calc_p0_p(self.mach, self.gamma)
        self.rho0_rho = IsentropicRelations.calc_rho0_rho(self.mach, self.gamma)
        self.a_aStar = IsentropicRelations.calc_A_Astar(self.mach, self.gamma)
        self.flowType = FS.SUPER_SONIC if self.mach > 1.0 else FS.SUB_SONIC

    def __str__(self) -> str:
        ff = "\nIsentropic Flow State at Mach: {}\n".format(round(self.mach, self._precision))
        ff2 = "---------------------------------------\n"
        first = "p0/p:  {}\n".format(round(self.p0_p, self._precision))
        second = "T0/T:  {}\n".format(round(self.t0_t, self._precision))
        third = "\u03C1_0/\u03C1: {}\n".format(round(self.rho0_rho, self._precision))
        fourth = "A/A*:  {}\n".format(round(self.a_aStar, self._precision))
        complete = "".join([ff, ff2, first, second, third, fourth])

        return "".join(
            [
                named_header("Isentropic Flow State at Mach", self.mach, precision=self._precision),
                seperator(),
                data_value_to_string("p0/p:", self.p0_p, precision=self._precision),
                data_value_to_string("T0/T:", self.t0_t, precision=self._precision, dot_line=True),
                data_value_to_string(
                    "{}0/{}:".format(lcg.rho, lcg.rho), self.rho0_rho, precision=self._precision
                ),
                data_value_to_string("A/A*:", self.a_aStar, precision=self._precision, dot_line=True),
                footer(),
            ]
        )

    def precision(self, precision: int) -> None:
        self._precision = int(precision)

    @staticmethod
    def calc_T0_T(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates T0_T for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        return 1 + (gamma - 1) / 2 * pow(mach, 2) - offset

    @staticmethod
    def calcMachFrom_T0_T(t0_t: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for T0_T """
        return sqrt((t0_t - 1) * 2 / (gamma - 1))

    @staticmethod
    def calc_p0_p(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates p0_p for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        return pow((1 + (gamma - 1) / 2 * pow(mach, 2)), gamma / (gamma - 1)) - offset

    @staticmethod
    def calcMachFrom_p0_p(p0_p: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for p0_p """

        return brenth(IsentropicRelations.calc_p0_p, 0, 30, args=(gamma, p0_p))

    @staticmethod
    def calc_rho0_rho(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates rho0_rho for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        return pow((1 + (gamma - 1) / 2 * pow(mach, 2)), 1 / (gamma - 1)) - offset

    @staticmethod
    def calcMachFrom_rho0_rho(rho0_rho: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for rho0_rho """

        return brenth(IsentropicRelations.calc_rho0_rho, 0, 30, args=(gamma, rho0_rho))

    @staticmethod
    def calc_A_Astar(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates A/A* for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        gm1 = gamma - 1
        gp1 = gamma + 1
        mSqr = pow(mach, 2)

        nonRaised = 2 / gp1 * (1 + gm1 / 2 * mSqr)

        return sqrt(pow(nonRaised, gp1 / gm1) / mSqr) - offset

    @staticmethod
    def calcMachFrom_A_Astar(A_Astar: float, gamma: float, flowType: FS = FS.SUPER_SONIC) -> float:
        """ Calcutes the given Mach associated for A/A* need to specify subsonic or supersonic solution """
        assert isinstance(flowType, FS)
        if A_Astar == 1.0:
            return A_Astar
        elif flowType == FS.SUPER_SONIC:
            return brenth(IsentropicRelations.calc_A_Astar, 1, 30, args=(gamma, A_Astar))
        elif flowType == FS.SUB_SONIC:
            return brenth(IsentropicRelations.calc_A_Astar, 0.001, 1, args=(gamma, A_Astar))
        else:
            print("Unsupported flow type passed to A/A* calculations. Type: {}".format(flowType.name))
            return nan
