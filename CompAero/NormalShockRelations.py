from enum import Enum
from math import sqrt, nan, pow, isnan
from _pytest.python_api import raises
from scipy.optimize import brenth

from CompAero.IsentropecRelations import IsentropicRelations
from CompAero.greek_letters import LowerCaseGreek as lcg
from CompAero.internal import (
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    checkValue,
    to_string,
    footer,
    named_header,
    seperator,
)

NORMAL_SHOCK_VALID_OPTIONS = [
    "gamma, mach",
    "gamma, P2/P1",
    "gamma, Rho2/Rho1",
    "gamma, T2/T1",
    "gamma, P02/P01",
    "gamma, P02/P1",
    "gamma, dwnStrm_mach",
]


class NORMAL_SHOCK_CHOICE(Enum):
    MACH = "gamma, mach"
    P2_P1 = "gamma, P2/P1"
    RHO2_RHO1 = "gamma, Rho2/Rho1"
    T2_T1 = "gamma, T2/T1"
    PO2_PO1 = "gamma, P02/P01"
    PO2_P1 = "gamma, P02/P1"
    M2 = "gamma, dwnStrm_mach"


class NormalShockRelations:
    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        p2_p1: float = nan,
        rho2_rho1: float = nan,
        t2_t1: float = nan,
        po2_po1: float = nan,
        po2_p1: float = nan,
        m2: float = nan,
    ) -> None:
        """ This class is a collective name space for basic calculations regarding Normal Shock Properties. 
        The constructor of this class can also determine the entire state of the flow given a partial state of the flow 

        Args:
            gamma (float): Ratio of specific heats
            mach (float, optional): mach number of the flow. Defaults to nan.
            p2_p1 (float, optional): Ratio of pressure behind shock wave to ahead of shock wave. Defaults to nan.
            rho2_rho1 (float, optional): Ratio of density behind shock wave to ahead of shock wave. Defaults to nan.
            t2_t1 (float, optional): Ratio of temperature behind shock wave to ahead of shock wave. Defaults to nan.
            po2_po1 (float, optional): Ratio of total pressure to behind shock wave to ahead of shock wave. Defaults to nan.
            po2_p1 (float, optional): Ratio Total pressure behind shock wave to static pressure ahead of shock wave. Defaults to nan.
            m2 (float, optional): mach number behind shock wave. Defaults to nan.
            
        Raises:
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and flow state cannot be determined
            
        Useage:
        To use this class pass gamma and one of the known parameters of the flow and the rest are calculated. 

        Valid_Combinations_of_Parameters:
            1: gamma, mach\n
            2: gamma, P2/P1\n
            3: gamma, Rho2/Rho1\n
            4: gamma, T2/T1\n
            5: gamma, P02/P01\n
            6: gamma, P02/P1\n
            7: gamma, dwnStrm_mach\n
        """
        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ Mach number of the flow """
        self.p2_p1 = p2_p1
        """ Ratio of pressure behind the shock wave to pressure before the shock wave P2/P1 """
        self.rho2_rho1 = rho2_rho1
        """ Ratio of density behind the shock wave to density before the shock wave rho2/rho1 """
        self.t2_t1 = t2_t1
        """ Ratio of temperature behind the shock wave to temperature before the shock wave T2/T1 """
        self.po2_po1 = po2_po1
        """ Ratio of pressure behind the shock wave to pressure before the shock wave P2/P1 """
        self.po2_p1 = po2_p1
        """ Ratio of total pressure behind the shock wave to pressure before the shock wave P02/P1 """
        self.mach2 = m2
        """ Mach number behind the shock wave """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        # Calculate parameters based on what was passed in
        if not checkValue(self.gamma):
            raise GammaNotDefinedError()

        if checkValue(self.mach):
            pass

        elif checkValue(self.p2_p1):
            self.mach = NormalShockRelations.calc_mach_from_p2_p1(self.p2_p1, self.gamma)

        elif checkValue(self.rho2_rho1):
            self.mach = NormalShockRelations.calc_mach_from_rho2_rho1(self.rho2_rho1, self.gamma)

        elif checkValue(self.t2_t1):
            self.mach = NormalShockRelations.calc_mach_from_T2_T1(self.t2_t1, self.gamma)

        elif checkValue(self.po2_po1):
            self.mach = NormalShockRelations.calc_mach_from_po2_po1(self.po2_po1, self.gamma)

        elif checkValue(self.po2_p1):
            self.mach = NormalShockRelations.calc_mach_from_po2_p1(self.po2_p1, self.gamma)

        elif checkValue(self.mach2):
            self.mach = NormalShockRelations.calc_mach_before_normal_shock_from_mach_after_shock(
                self.mach2, self.gamma
            )

        else:
            raise InvalidOptionCombinationError()

        if checkValue(self.mach):
            if self.mach >= 1.0:
                self.__calculateStateFromMach()
            else:
                self.mach = nan

    def __calculateStateFromMach(self) -> None:
        self.p2_p1 = NormalShockRelations.calc_p2_p1(self.mach, self.gamma)
        self.rho2_rho1 = NormalShockRelations.calc_rho2_rho1(self.mach, self.gamma)
        self.t2_t1 = NormalShockRelations.calc_T2_T1(self.mach, self.gamma)
        self.po2_po1 = NormalShockRelations.calc_po2_po1(self.mach, self.gamma)
        self.po2_p1 = NormalShockRelations.calc_po2_p1(self.mach, self.gamma)
        self.mach2 = NormalShockRelations.calc_mach_after_normal_shock(self.mach, self.gamma)

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Normal Shock Relations at Mach", self.mach, self.precision),
                seperator(),
                to_string("P2/P1", self.p2_p1, self.precision),
                to_string("{}2/{}1".format(*[lcg.rho] * 2), self.rho2_rho1, self.precision, dot_line=True),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("P02/P1", self.po2_p1, self.precision),
                to_string("Dowstream Mach", self.mach2, self.precision, dot_line=True),
                footer(),
            ]
        )

    @staticmethod
    def calc_p2_p1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of Pressure behind the shock wave to that of before the shock wave P2/P1

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Raises:
            ValueError: Raised if mach < 1.0
            
        Returns:
            float: P2/P1
        """

        if mach < 1.0:
            raise ValueError("Normal Shocks Require a mach greater than 1")

        return 1 + 2 * gamma / (gamma + 1) * (pow(mach, 2) - 1) - offset

    @staticmethod
    def calc_mach_from_p2_p1(p2_p1: float, gamma: float) -> float:
        """Calculates the mach number based on the ratio Pressure behind the shock wave to that of before the shock wave P2/P1

        Args:
            p2_p1 (float): Ratio of Pressure behind the shock wave to that of before the shock wave P2/P1
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(NormalShockRelations.calc_p2_p1, 1.0, 30, args=(gamma, p2_p1))

    @staticmethod
    def calc_rho2_rho1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of density behind the shock wave to that of density before the shock wave Rho2/Rho1

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Raises:
            ValueError: Raised if mach < 1.0
            
        Returns:
            float: Rho2/Rho1
        """
        if mach < 1.0:
            raise ValueError("Normal Shocks Require a mach greater than 1")

        return (gamma + 1) * pow(mach, 2) / (2 + (gamma - 1) * pow(mach, 2)) - offset

    @staticmethod
    def calc_mach_from_rho2_rho1(rho2_rho1: float, gamma: float) -> float:
        """Calculates the mach number based on the ratio density behind the shock wave to that of before the shock wave Rho2/Rho1

        Args:
            p2_p1 (float): Ratio of density behind the shock wave to that of before the shock wave Rho2/Rho1
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(NormalShockRelations.calc_rho2_rho1, 1.0, 30, args=(gamma, rho2_rho1))

    @staticmethod
    def calc_T2_T1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of Temperature behind the shock wave to that of Temperature before the shock wave T2/T1

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Raises:
            ValueError: Raised if mach < 1.0
            
        Returns:
            float: T2/T1
        """
        if mach < 1.0:
            raise ValueError("Normal Shocks Require a mach greater than 1")

        return (
            NormalShockRelations.calc_p2_p1(mach, gamma) / NormalShockRelations.calc_rho2_rho1(mach, gamma)
            - offset
        )

    @staticmethod
    def calc_mach_from_T2_T1(t2_t1: float, gamma: float) -> float:
        """Calculates the mach number based on the ratio of temperature behind the shock wave to that of before the shock wave T2/T1

        Args:
            t2_t1 (float): Ratio of temperature behind the shock wave to that of before the shock wave T2/T1
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(NormalShockRelations.calc_T2_T1, 1, 30, args=(gamma, t2_t1))

    @staticmethod
    def calc_mach_after_normal_shock(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates the mach number behind the shock wave

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Raises:
            ValueError: Raised if mach < 1.0
            
        Returns:
            float: mach number behind shock wave
        """
        if mach < 1.0:
            raise ValueError("Normal Shocks Require a mach greater than 1")
        gm1 = gamma - 1
        mSqr = pow(mach, 2)
        num = 1 + gm1 / 2 * mSqr
        denom = gamma * mSqr - gm1 / 2
        return sqrt(num / denom) - offset

    @staticmethod
    def calc_mach_before_normal_shock_from_mach_after_shock(mach2: float, gamma: float) -> float:
        """Calculates the mach number before a shock wave based on the mach number after a shock wave

        Args:
            mach2 (float): mach number after the shock wave
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(NormalShockRelations.calc_mach_after_normal_shock, 1, 30, args=(gamma, mach2))

    @staticmethod
    def calc_po2_po1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of Total Pressure behind the shock wave to that of Total Pressure before the shock wave P02/P01

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Raises:
            ValueError: Raised if mach < 1.0
            
        Returns:
            float: P02/P01
        """
        if mach < 1.0:
            raise ValueError("Normal Shocks Require a mach greater than 1")
        m2 = NormalShockRelations.calc_mach_after_normal_shock(mach, gamma)
        upstream = IsentropicRelations(gamma, mach=mach)
        downstream = IsentropicRelations(gamma, mach=m2)
        p2_p1 = NormalShockRelations.calc_p2_p1(mach, gamma)
        return downstream.p0_p * p2_p1 / upstream.p0_p - offset

    @staticmethod
    def calc_mach_from_po2_po1(po2_po1: float, gamma: float) -> float:
        """Calculates the mach number based on the ratio of Total Pressure behind the shock wave to that of before the shock wave P02/P01

        Args:
            po2_po1 (float): Ratio of Total Pressure behind the shock wave to that of before the shock wave P02/P01
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(NormalShockRelations.calc_po2_po1, 1, 30, args=(gamma, po2_po1))

    @staticmethod
    def calc_po2_p1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculates Ratio of Total Pressure behind the shock wave to that of static Pressure before the shock wave P02/P1

        Args:
            mach (float): mach number of the flow
            gamma (float): ratio of specific heats
            offset (float, optional): offset that can be used for root finding for a specific value. Defaults to 0.0.

        Raises:
            ValueError: Raised if mach < 1.0
            
        Returns:
            float: P02/P1
        """
        gm1 = gamma - 1
        gp1 = gamma + 1
        mSqr = pow(mach, 2)

        firstratio = pow(gp1, 2) * mSqr / (4 * gamma * mSqr - 2 * gm1)
        secondRatio = (1 - gamma + 2 * gamma * mSqr) / gp1
        return pow(firstratio, gamma / gm1) * secondRatio - offset

    @staticmethod
    def calc_mach_from_po2_p1(po2_p1: float, gamma: float) -> float:
        """Calculates the mach number based on the ratio of Total Pressure behind the shock wave to that of before the shock wave P02/P1
        Also known as the "Rayleigh Pitot Tube Formula"

        Args:
            po2_p1 (float): Ratio of temperature behind the shock wave to that of before the shock wave P02/P1
            gamma (float): ratio of specific heats

        Returns:
            float: mach number
        """
        return brenth(NormalShockRelations.calc_po2_p1, 1, 30, args=(gamma, po2_p1))
