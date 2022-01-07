from math import atan, sqrt, nan, pow, radians, degrees
from types import DynamicClassAttribute
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from CompAero.internal import checkValue, footer, named_header, named_subheader, seperator, to_string
from CompAero.ObliqueShockRelations import ObliqueShockRelations
from CompAero.greek_letters import LowerCaseGreek as lcg, Misc


class PrandtlMeyer:
    def __init__(
        self,
        gamma: float,
        mach: float = nan,
        nu: float = nan,
        mu: float = nan,
        dwnstreamNu: float = nan,
        dwnStreamMu: float = nan,
        dwnStreamMach: float = nan,
        deflectionAngle: float = nan,
        inDegrees: bool = True,
    ) -> None:
        self.__gamma = gamma
        self.__mach = mach
        self.__nu = nu
        self.__mu = mu
        self.__deflectionAngle = deflectionAngle
        self.__degrees = inDegrees
        self.__preciscion = 4

        self.__dwnStrmNu = dwnstreamNu
        self.__dwnStrmMu = dwnStreamMu
        self.__dwnStrmMach = dwnStreamMach

        if not self.__degrees:
            self.__deflectionAngle = degrees(self.__deflectionAngle)

        if not checkValue(self.__gamma):
            raise ValueError(
                "Valid gamma is required for Prandtl Meyer calculations. Given {}".format(self.__gamma)
            )

        if checkValue(self.__mach):
            pass
        elif checkValue(self.__nu):
            self.__mach = PrandtlMeyer.calcMachGivenNu(self.__nu, self.__gamma)

        elif checkValue(self.__mu):
            self.__mach = ObliqueShockRelations.calcMachFromMachWaveAngle(radians(self.__mu))

        elif checkValue(self.__deflectionAngle) and checkValue(self.__dwnStrmMach):
            self.__dwnStrmNu = PrandtlMeyer.calcNuFromMach(self.__dwnStrmMach, self.__gamma)
            self.__nu = self.__dwnStrmNu - self.__deflectionAngle
            self.__mach = PrandtlMeyer.calcMachGivenNu(self.__nu, self.__gamma)

        elif checkValue(self.__deflectionAngle) and checkValue(self.__dwnStrmNu):
            self.__nu = self.__dwnStrmNu - self.__deflectionAngle
            self.__mach = PrandtlMeyer.calcMachGivenNu(self.__nu, self.__gamma)

        elif checkValue(self.__deflectionAngle) and checkValue(self.__dwnStrmMu):
            self.__dwnStrmMach = ObliqueShockRelations.calcMachFromMachWaveAngle(radians(self.__dwnStrmMu))
            self.__dwnStrmNu = PrandtlMeyer.calcNuFromMach(self.__dwnStrmMach, self.__gamma)
            self.__nu = self.__dwnStrmNu - self.__deflectionAngle
            self.__mach = PrandtlMeyer.calcMachGivenNu(self.__nu, self.__gamma)

        else:
            raise ValueError(
                Fore.BLACK
                + Back.RED
                + "Incorrect Combination of input arguments given to {}".format(self.__class__.__name__)
                + Style.RESET_ALL
                + "\t"
            )

        if checkValue(self.__mach) and self.__mach >= 1.0:
            self.__calculateState()

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Prandtl Relations at Mach", self.mach, self.precision),
                seperator(),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string(lcg.nu, self.nu, self.precision, dot_line=True),
                to_string(lcg.mu, self.mu, self.precision),
                seperator(),
                named_subheader("Downstream Conditions"),
                to_string("Mach", self.__dwnStrmMach, self.precision),
                to_string(lcg.nu, self.__dwnStrmNu, self.precision, dot_line=True),
                to_string(lcg.mu, self.downStreamMu, self.precision),
                to_string(
                    "Flow Deflection Angle [{}]".format(lcg.theta),
                    self.__deflectionAngle,
                    self.precision,
                    dot_line=True,
                ),
                footer(),
            ]
        )

    def __calculateState(self) -> None:
        self.__nu = PrandtlMeyer.calcNuFromMach(self.__mach, self.__gamma)
        self.__mu = degrees(ObliqueShockRelations.calcMachWaveAngle(self.__mach))

        if not checkValue(self.__deflectionAngle):
            return

        self.__dwnStrmNu = self.__deflectionAngle + self.__nu
        self.__dwnStrmMach = PrandtlMeyer.calcMachGivenNu(self.__dwnStrmNu, self.__gamma)
        self.__dwnStrmMu = degrees(ObliqueShockRelations.calcMachWaveAngle(self.__dwnStrmMach))

    @property
    def gamma(self) -> float:
        """Ratio of specific heats"""
        return self.__gamma

    @property
    def mach(self) -> float:
        """ Free stream mach number"""
        return self.__mach

    @property
    def nu(self) -> float:
        """Prandtl Meyer Result from free stream mach number"""
        return self.__nu

    @property
    def mu(self) -> float:
        """ Mach wave angle associated with free stream mach number"""
        return self.__mu

    @property
    def deflectionAngle(self) -> float:
        """ Angle the flow is deflected by """
        return self.__deflectionAngle

    @property
    def downStreamNu(self) -> float:
        """ Returns the result of the prandtl meyer function of the downstream mach number"""
        return self.__dwnStrmNu

    @property
    def downStreamMu(self) -> float:
        """ Returns the mach wave angle associated with the down stream mach number"""
        return self.__dwnStrmMu

    @property
    def downStreamMach(self) -> float:
        """Returns the down stream mach number"""
        return self.__dwnStrmMach

    @property
    def precision(self) -> int:
        """ Precision of any printed strings"""
        return self.__preciscion

    @precision.setter
    def precision(self, precision) -> None:
        """Sets the preciscion for any printed strings"""
        self.__preciscion = precision

    @staticmethod
    def calcNuFromMach(mach: float, gamma: float, offset: float = 0.0) -> float:
        """Calculated Nu given a mach number. Offset can be used for root finding"""
        if mach <= 1.0:
            return 0.0

        gp1 = gamma + 1
        gm1 = gamma - 1
        mSqrMinus1 = pow(mach, 2) - 1
        return degrees(sqrt(gp1 / gm1) * atan(sqrt(gm1 / gp1 * mSqrMinus1)) - atan(sqrt(mSqrMinus1))) - offset

    @staticmethod
    def calcMachGivenNu(nu: float, gamma: float) -> float:
        """ Calculates Mach number for a given Nu"""
        if nu <= 0.0:
            return 1.0

        return brenth(PrandtlMeyer.calcNuFromMach, 1 + 1e-9, 30, args=(gamma, nu))

