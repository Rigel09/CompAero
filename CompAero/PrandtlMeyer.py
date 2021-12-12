from math import atan, sqrt, nan, pow, radians, degrees
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from CompAero.common import checkValue
from CompAero.ObliqueShockRelations import ObliqueShockRelations
from GreekLetters.greekLetters import LowerCaseGreek as lcg, Misc


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

        if checkValue(self.__mach):
            self.__calculateState()

    def __str__(self) -> str:
        gammaStr = str(round(self.__gamma, self.__preciscion))
        machStr = str(round(self.__mach, self.__preciscion))
        nuStr = str(round(self.__nu, self.__preciscion))
        muStr = str(round(self.__mu, self.__preciscion))
        dwnNuStr = str(round(self.__dwnStrmNu, self.__preciscion))
        dwnMuStr = str(round(self.__dwnStrmMu, self.__preciscion))
        deflectionStr = str(round(self.__deflectionAngle, self.__preciscion))
        dwnMachStr = str(round(self.__dwnStrmMach, self.__preciscion))

        width = 50 - 2  # Width minus edges
        sep = "|{:{width}}|\n".format("", width=width)

        # Initial Conditions
        header = "|{:=^{width}}|\n".format("Prandtl Relations at Mach: {}".format(machStr), width=width)
        gammaS = "|{:<{width}}{}|\n".format(lcg.gamma, gammaStr, width=width - len(gammaStr))
        nuS = "|{:-<{width}}{}|\n".format(lcg.nu, nuStr, width=width - len(nuStr))
        muS = "|{:<{width}}{}|\n".format(lcg.mu, muStr, width=width - len(muStr))

        jumpSep = "|{:-^{width}}|\n".format("Downstream Conditions", width=width)
        dwnMachS = "|{:-<{width}}{}|\n".format("Mach", dwnMachStr, width=width - len(dwnMachStr))
        dwnNu = "|{:<{width}}{}|\n".format(lcg.nu, dwnNuStr, width=width - len(dwnNuStr))
        dwnMu = "|{:-<{width}}{}|\n".format(lcg.mu, dwnMuStr, width=width - len(dwnMuStr))
        dwnDeflection = "|{:<{width}}{}|\n".format(
            "Flow Deflection Angle [" + lcg.theta + Misc.degreeSym + "]",
            deflectionStr,
            width=width - len(deflectionStr),
        )

        finish = "|{:=^{width}}|\n".format("", width=width)

        return "".join(
            [header, sep, gammaS, nuS, muS, sep, jumpSep, dwnMachS, dwnNu, dwnMu, dwnDeflection, finish]
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
        gp1 = gamma + 1
        gm1 = gamma - 1
        mSqrMinus1 = pow(mach, 2) - 1
        return degrees(sqrt(gp1 / gm1) * atan(sqrt(gm1 / gp1 * mSqrMinus1)) - atan(sqrt(mSqrMinus1))) - offset

    @staticmethod
    def calcMachGivenNu(nu: float, gamma: float) -> float:
        """ Calculates Mach number for a given Nu"""
        return brenth(PrandtlMeyer.calcNuFromMach, 1 + 1e-9, 30, args=(gamma, nu))

