from math import sqrt, nan, pow, isnan
from scipy.optimize import brenth

from CompAero.IsentropecRelations import IsentropicRelations


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
        self.gamma = gamma
        self.mach = mach
        self.p2_p1 = p2_p1
        self.rho2_rho1 = rho2_rho1
        self.t2_t1 = t2_t1
        self.po2_po1 = po2_po1
        self.po2_p1 = po2_p1
        self.mach2 = m2
        self._precision = 4

        # Calculate parameters based on what was passed in
        if isnan(self.gamma) or self.gamma < 0.0:
            return

        if self.__checkValue(self.p2_p1):
            self.mach = NormalShockRelations.calcMachFrom_p2_p1(self.p2_p1, self.gamma)

        elif self.__checkValue(self.rho2_rho1):
            self.mach = NormalShockRelations.calcMachFrom_rho2_rho1(self.rho2_rho1, self.gamma)

        elif self.__checkValue(self.t2_t1):
            self.mach = NormalShockRelations.calcMachFrom_T2_T1(self.t2_t1, self.gamma)

        elif self.__checkValue(self.po2_po1):
            self.mach = NormalShockRelations.calcMachFrom_po2_po1(self.po2_po1, self.gamma)

        elif self.__checkValue(self.po2_p1):
            self.mach = NormalShockRelations.calcMachFrom_po2_p1(self.po2_p1, self.gamma)

        elif self.__checkValue(self.mach2):
            self.mach = NormalShockRelations.calcMachFrom_mach2(self.mach2, self.gamma)

        if self.__checkValue(self.mach):
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
        self.mach2 = NormalShockRelations.calc__mach2(self.mach, self.gamma)

    def __checkValue(self, var: float) -> bool:
        if isnan(var):
            return False

        if var < 0:
            return False

        return True

    def setPrecision(self, precision: int) -> None:
        self._precision = int(precision)

    def __str__(self) -> str:
        mach1Str = str(round(self.mach, self._precision))
        pRatioStr = str(round(self.p2_p1, self._precision))
        rhoRatioStr = str(round(self.rho2_rho1, self._precision))
        tempRatioStr = str(round(self.t2_t1, self._precision))
        poRatioStr = str(round(self.po2_po1, self._precision))
        poprRatioStr = str(round(self.po2_p1, self._precision))
        mach2Str = str(round(self.mach2, self._precision))

        width = 50 - 2  # Width minus edges

        header = "|{:=^{width}}|\n".format("Normal Shock Relations at Mach: {}".format(mach1Str), width=width)
        pRatio = "|{:<{width}}{}|\n".format("p2/p1:", pRatioStr, width=width - len(pRatioStr))
        rhoRatio = "|{:-<{width}}{}|\n".format(
            "\u03C12/\u03C11:", rhoRatioStr, width=width - len(rhoRatioStr)
        )
        tempRatio = "|{:<{width}}{}|\n".format("T2/T1:", tempRatioStr, width=width - len(tempRatioStr))
        poRatio = "|{:-<{width}}{}|\n".format("po2/po1:", poRatioStr, width=width - len(poRatioStr))
        poprRatio = "|{:<{width}}{}|\n".format("po2/p1:", poprRatioStr, width=width - len(poprRatioStr))
        mach2 = "|{:<{width}}{}|\n".format("Downstream Mach:", mach2Str, width=width - len(mach2Str))
        sep = "|{:{width}}|\n".format("", width=width)
        finish = "|{:=^{width}}|".format("", width=width)

        return "".join([header, sep, pRatio, rhoRatio, tempRatio, poRatio, poprRatio, mach2, finish,])

    @staticmethod
    def calc_p2_p1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates p2/p1 for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        return 1 + 2 * gamma / (gamma + 1) * (pow(mach, 2) - 1) - offset

    @staticmethod
    def calcMachFrom_p2_p1(p2_p1: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for p2_p1 """
        return brenth(NormalShockRelations.calc_p2_p1, 0.001, 30, args=(gamma, p2_p1))

    @staticmethod
    def calc_rho2_rho1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates p2/p1 for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        return (gamma + 1) * pow(mach, 2) / (2 + (gamma - 1) * pow(mach, 2)) - offset

    @staticmethod
    def calcMachFrom_rho2_rho1(rho2_rho1: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for rho2/rho1 """
        return brenth(NormalShockRelations.calc_rho2_rho1, 0.001, 30, args=(gamma, rho2_rho1))

    @staticmethod
    def calc_T2_T1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates T2/T1 for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """

        return (
            NormalShockRelations.calc_p2_p1(mach, gamma) / NormalShockRelations.calc_rho2_rho1(mach, gamma)
            - offset
        )

    @staticmethod
    def calcMachFrom_T2_T1(t2_t1: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for t2/t1 """
        return brenth(NormalShockRelations.calc_T2_T1, 1, 30, args=(gamma, t2_t1))

    @staticmethod
    def calc__mach2(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates downstream mach for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """
        gm1 = gamma - 1
        mSqr = pow(mach, 2)
        num = 1 + gm1 / 2 * mSqr
        denom = gamma * mSqr - gm1 / 2
        return sqrt(num / denom) - offset

    @staticmethod
    def calcMachFrom_mach2(mach2: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for downstream mach """
        return brenth(NormalShockRelations.calc__mach2, 1, 30, args=(gamma, mach2))

    @staticmethod
    def calc_po2_po1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates po2/po1 for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
        """
        m2 = NormalShockRelations.calc__mach2(mach, gamma)
        upstream = IsentropicRelations(gamma, mach=mach)
        downstream = IsentropicRelations(gamma, mach=m2)
        p2_p1 = NormalShockRelations.calc_p2_p1(mach, gamma)
        return downstream.p0_p * p2_p1 / upstream.p0_p - offset

    @staticmethod
    def calcMachFrom_po2_po1(po2_po1: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for po2/po1 """
        return brenth(NormalShockRelations.calc_po2_po1, 1, 30, args=(gamma, po2_po1))

    @staticmethod
    def calc_po2_p1(mach: float, gamma: float, offset: float = 0.0) -> float:
        """ 
            Calculates po2/p1 for a give Mach and gamma combination
            Can be used for root finding if a given off set is applied
            [Rayleigh Pitot Tube Formula] 
        """
        gm1 = gamma - 1
        gp1 = gamma + 1
        mSqr = pow(mach, 2)

        firstratio = pow(gp1, 2) * mSqr / (4 * gamma * mSqr - 2 * gm1)
        secondRatio = (1 - gamma + 2 * gamma * mSqr) / gp1
        return pow(firstratio, gamma / gm1) * secondRatio - offset

    @staticmethod
    def calcMachFrom_po2_p1(po2_p1: float, gamma: float) -> float:
        """ Calcutes the given Mach associated for po2/p1  [Rayleigh Pitot Tube Formula] """
        return brenth(NormalShockRelations.calc_po2_p1, 1, 30, args=(gamma, po2_p1))
