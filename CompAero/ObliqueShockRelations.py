from math import sqrt, nan, pow, radians, sin, cos, atan, tan, asin, degrees, asin
from scipy.optimize import brenth
import matplotlib.pyplot as plt
import numpy as np

from CompAero.NormalShockRelations import NormalShockRelations
from CompAero.internal import (
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    checkValue,
    footer,
    named_header,
    named_subheader,
    seperator,
    to_string,
    ShockType,
)
from CompAero.greek_letters import LowerCaseGreek as lcg

# TODO: Subclassing from the normal shock relations doesnt seem to make as much sense as it used to
# get rid of the subclassing and mach ObliqueShockRelations it's own class
class ObliqueShockRelations(NormalShockRelations):
    """ This class is a collective name space for basic calculations regarding Oblique Shock Properties. 
    The constructor of this class can also determine the entire state of the flow given a partial state of the flow 

    Args:
        gamma (float): Ratio of specific heats
        useDegrees (bool, optional): Angles provided to class are in degrees. Output will still be in radians. Defaults to True.
        shockAngle (float, optional): Angle of the oblique shock wave. Defaults to nan.
        wedgeAngle (float, optional): Angle of the wedge that is deflecting the flow. Defaults to nan.
        mn1 (float, optional): Normal component of mach number ahead of the shock. Defaults to nan.
        mn2 (float, optional): Normal component of mach number behind the shock. Defaults to nan.
        mach (float, optional): Mach number ahead of the shock. Defaults to nan.
        p2_p1 (float, optional): Ratio of pressure behind shock wave to ahead of shock wave. Defaults to nan.
        rho2_rho1 (float, optional): Ratio of density behind shock wave to ahead of shock wave. Defaults to nan.
        t2_t1 (float, optional): Ratio of temperature behind shock wave to ahead of shock wave. Defaults to nan.
        po2_po1 (float, optional): Ratio of total pressure to behind shock wave to ahead of shock wave. Defaults to nan.
        po2_p1 (float, optional): Ratio Total pressure behind shock wave to static pressure ahead of shock wave. Defaults to nan.
        m2 (float, optional): mach number behind shock wave. Defaults to nan.
        shockType (ShockType, optional): Describes whether the oblique shock is a Week or Strong shock. Defaults to ShockType.WEAK.

    Raises:
        RuntimeError: If shock wave angle is found to be less than the mach wave angle at
        RuntimeError: If the deflection angle is found to be greater than the max deflection angle for the flow
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and flow state cannot be determined
        
    Valid_Combinations_of_Parameters:
        1: gamma, mach\n
        2: gamma, mach normal 1, shock angle \n
        3: gamma, mach behind shock, wedge angle, shock angle\n
        4: gamma, shock angle, wedge angle\n
        5: gamma, P2/P1, shock angle\n
        6: gamma, Rho2/Rho1, shock angle\n
        7: gamma, T2/T1, shock angle\n
        8: gamma, P02/P01, shock angle\n
        9: gamma, P02/P1, shock angle\n
        10: gamma, dwnStrm_mach, shock angle\n
    """

    def __init__(
        self,
        gamma: float,
        useDegrees: bool = True,
        shockAngle: float = nan,
        wedgeAngle: float = nan,
        mn1: float = nan,
        mn2: float = nan,
        mach: float = nan,
        p2_p1: float = nan,
        rho2_rho1: float = nan,
        t2_t1: float = nan,
        po2_po1: float = nan,
        po2_p1: float = nan,
        m2: float = nan,
        shockType: ShockType = ShockType.WEAK,
    ) -> None:
        self.useDegrees = useDegrees
        """ Angles passed to the constructor are in degrees"""
        self.shockAngle = shockAngle
        """ Angle Oblique shock"""
        self.wedgeAngle = wedgeAngle
        """ Angle of the flow deflection (wedge) angle"""
        self.machNorm1 = mn1
        """ The normal component of the mach number ahead the shock wave"""
        self.machNorm2 = mn2
        """ The normal component of the mach number behind the shock wave"""
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
        self.shockType: ShockType = shockType
        """ The type of shock whether it is strong or weak """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        if self.useDegrees and checkValue(self.shockAngle):
            self.shockAngle = radians(self.shockAngle)
        if self.useDegrees and checkValue(self.wedgeAngle):
            self.wedgeAngle = radians(self.wedgeAngle)

        if not checkValue(self.gamma):
            raise GammaNotDefinedError()

        if checkValue(self.mach):
            pass

        elif checkValue(self.machNorm1) and checkValue(self.shockAngle):
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue(m2) and checkValue(self.wedgeAngle) and checkValue(self.shockAngle):
            self.machNorm2 = ObliqueShockRelations.calc_mach_normal_behind_shock_from_mach_behind_shock(
                m2, self.shockAngle, self.wedgeAngle
            )
            if self.machNorm2 > 1:
                assert ValueError("Normal component of downstream mach number has to be less than 1")
            self.machNorm1 = NormalShockRelations.calc_mach_before_normal_shock_from_mach_after_shock(
                self.machNorm2, self.gamma
            )
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue(self.shockAngle) and checkValue(self.wedgeAngle):
            self.mach = ObliqueShockRelations.calc_mach_from_theta_beta_mach(
                self.shockAngle, self.wedgeAngle, self.gamma
            )
            machAngle = ObliqueShockRelations.calc_mach_wave_angle(self.mach)
            maxShockAngle = ObliqueShockRelations.calc_max_shock_angle(self.mach, self.gamma)
            maxDeflectionAngle = ObliqueShockRelations.calc_max_flow_deflection_angle(
                maxShockAngle, self.mach, self.gamma
            )

            if machAngle > self.shockAngle:
                raise RuntimeError(
                    "Shock Angle, {}, cannot be less than that of Mach Wave, {}".format(
                        round(self.shockAngle, 4), round(machAngle, 4)
                    )
                )

            if maxDeflectionAngle < self.wedgeAngle:
                raise RuntimeError(
                    "SHOCK DETACHED!! Deflection Angles: Given: {}\tMax: {}".format(
                        round(self.wedgeAngle, 4), round(maxDeflectionAngle, 4)
                    )
                )

        elif checkValue([po2_p1, self.shockAngle]):
            self.machNorm1 = NormalShockRelations.calc_mach_from_po2_p1(po2_p1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue([p2_p1, self.shockAngle]):
            self.machNorm1 = NormalShockRelations.calc_mach_from_p2_p1(p2_p1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue([rho2_rho1, self.shockAngle]):
            self.machNorm1 = NormalShockRelations.calc_mach_from_rho2_rho1(rho2_rho1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue([t2_t1, self.shockAngle]):
            self.machNorm1 = NormalShockRelations.calc_mach_from_T2_T1(t2_t1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue([po2_po1, self.shockAngle]):
            self.machNorm1 = NormalShockRelations.calc_mach_from_po2_po1(po2_po1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        elif checkValue(self.machNorm2) and checkValue(self.shockAngle):
            self.machNorm1 = NormalShockRelations.calc_mach_before_normal_shock_from_mach_after_shock(
                self.machNorm2, self.gamma
            )
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.machNorm1, self.shockAngle
            )

        else:
            raise InvalidOptionCombinationError()

        if checkValue(self.mach) and ((checkValue(self.shockAngle) or checkValue(self.wedgeAngle))):
            self.__calculateState()

        super().__init__(self.gamma, mach=self.machNorm1)
        self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
            self.machNorm1, self.shockAngle
        )
        self.machNorm2 = self.mach2
        self.mach2 = ObliqueShockRelations.calc_mach_behind_shock(
            self.machNorm2, self.wedgeAngle, self.shockAngle
        )

    def __calculateState(self) -> None:
        if checkValue(self.wedgeAngle):
            if self.shockType == ShockType.WEAK:
                self.shockAngle = ObliqueShockRelations.calc_beta_from_theta_beta_mach_weak(
                    self.wedgeAngle, self.mach, self.gamma
                )

            elif self.shockType == ShockType.STRONG:
                self.shockAngle = ObliqueShockRelations.calc_beta_from_theta_beta_mach_strong(
                    self.wedgeAngle, self.mach, self.gamma
                )

            else:
                raise ValueError(
                    "Incorrect shock type specified -> [{}]. Choices: Strong, Weak".format(
                        self.shockType.name
                    )
                )

        elif checkValue(self.shockAngle):
            self.wedgeAngle = ObliqueShockRelations.calc_theta_from_theta_beta_mach(
                self.shockAngle, self.mach, self.gamma
            )

        self.machNorm1 = ObliqueShockRelations.calc_mach_normal_ahead_shock(self.mach, self.shockAngle)

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Oblique Shock Relations at Mach", self.mach, self.precision),
                seperator(),
                named_subheader("Upstream Conditions"),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string("Mach", self.mach, self.precision, dot_line=True),
                to_string("Mach Normal Component", self.machNorm1, self.precision),
                to_string(
                    "Flow Deflection Angle {}".format(lcg.theta),
                    self.wedgeAngle,
                    self.precision,
                    dot_line=True,
                ),
                to_string("Shock Angle {}".format(lcg.beta), self.shockAngle, self.precision),
                to_string("Flow Turn Type", self.shockType.name, self.precision),
                seperator(),
                named_subheader("Shock Jump Conditions"),
                to_string("P2/P1", self.p2_p1, self.precision),
                to_string("{}2/{}1".format(*[lcg.rho] * 2), self.rho2_rho1, self.precision, dot_line=True),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("P02/P1", self.po2_p1, self.precision),
                seperator(),
                named_subheader("Downstream Conditions"),
                to_string("Mach", self.mach2, self.precision, dot_line=True),
                to_string("Mach Normal Component", self.machNorm2, self.precision),
                footer(),
            ]
        )

    @staticmethod
    def calc_mach_normal_ahead_shock(mach: float, beta: float) -> float:
        """ Calculates the normal component of the mach number ahead of the shock wave

        Args:
            mach (float): Mach number of flow ahead of shock wave
            beta (float): Angle of oblique shock in radians

        Raises:
            ValueError: Raised if mach number is less than 1.0

        Returns:
            float: Normal component of upstream mach number
        """
        if mach < 1.0:
            raise ValueError("Normal Shocks Require a mach greater than 1")

        machWave = ObliqueShockRelations.calc_mach_wave_angle(mach)

        if abs(beta - machWave) < 1e-5:
            return 1.0

        if beta < machWave:
            return nan

        return mach * sin(beta)

    @staticmethod
    def calc_mach_ahead_shock_from_mach_normal_ahead_shock(machNormal1: float, beta: float) -> float:
        """ Calculates the upstream mach number from the normal component of the upstream mach number

        Args:
            machNormal1 (float): Normal Component of mach number ahead of the shock wave
            beta (float): Angle of oblique shock in radians

        Returns:
            float: Returns value of mach number of the flow ahead of the shock wave
        """
        return machNormal1 / sin(beta)

    @staticmethod
    def calc_beta_from_mach_mach_normal_ahead_shock(mach: float, machNormal1: float) -> float:
        """ Calculates the Oblique shock angle from the normal component of the mach number that is ahead of the shock wave

        Args:
            mach (float): Mach number of the flow ahead of the shock wave
            machNormal1 (float): Normal Component of mach number of the flow that is ahead of the oblique shock wave

        Returns:
            float: Oblique Shock Angle
        """
        return asin(machNormal1 / mach)

    @staticmethod
    def calc_mach_behind_shock(machNormal2: float, theta: float, beta: float) -> float:
        """ Calculates the Mach number behind the oblique shock wave

        Args:
            machNormal2 (float): Normal Component of the mach number behind the shock wave
            theta (float): Flow Deflection Angle (radians) (Wedge angle)
            beta (float): Oblique shock angle (radians)

        Returns:
            float: Mach number of flow behind the oblique shock
        """
        return machNormal2 / sin(beta - theta)

    @staticmethod
    def calc_mach_normal_behind_shock_from_mach_behind_shock(
        mach2: float, beta: float, theta: float
    ) -> float:
        """ Calculates the normal component of the mach number behind the oblique shock 

        Args:
            mach2 (float): Mach number of flow behind the oblique shock
            beta (float): Oblique shock angle (radians)
            theta (float): Flow deflections (Wedge) angle (radians)

        Returns:
            float: Normal component of mach numnber of the flow behind the shock wave
        """
        return mach2 * sin(beta - theta)

    @staticmethod
    def calc_theta_from_theta_beta_mach(beta: float, mach: float, gamma: float, offset: float = 0.0) -> float:
        """ Impliments the Theta-Beta-Mach (TBM) equation. Solves for Theta

        Args:
            beta (float): Oblique shock angle (radians)
            mach (float): Mach number of flow ahead of the shock wave
            gamma (float): Ratio of specific heats
            offset (float, optional): [description]. Defaults to 0.0.

        Returns:
            float: Flow deflection (Wedge) angle (radians)
        """
        mSqr = pow(mach, 2)
        num = mSqr * pow(sin(beta), 2) - 1
        denom = mSqr * (gamma + cos(2 * beta)) + 2
        theta = atan(2 * 1 / tan(beta) * num / denom) - offset
        return theta

    @staticmethod
    def calc_beta_from_theta_beta_mach_weak(theta: float, mach: float, gamma: float) -> float:
        """ Impliments the Theta-Beta-Mach (TBM) equation. Solves for Beta (shock angle) assuming the shock is weak

        Args:
            theta (float): Flow deflection (Wedge) angle (radians)
            mach (float): Mach number of the flow ahead of the oblique shock
            gamma (float): ratio of specific heats

        Returns:
            float: Oblique shock angle (radians)
        """
        maxShockAngle = ObliqueShockRelations.calc_max_shock_angle(mach, gamma)
        minShockAngle = ObliqueShockRelations.calc_mach_wave_angle(mach)
        return brenth(
            ObliqueShockRelations.calc_theta_from_theta_beta_mach,
            minShockAngle,
            maxShockAngle,
            args=(mach, gamma, theta),
        )

    @staticmethod
    def calc_beta_from_theta_beta_mach_strong(theta: float, mach: float, gamma: float) -> float:
        """ Impliments the Theta-Beta-Mach (TBM) equation. Solves for Beta (shock angle) assuming a strong shock wave

        Args:
            theta (float): Flow deflection (Wedge) angle (radians)
            mach (float): Mach number of the flow ahead of the oblique shock
            gamma (float): ratio of specific heats

        Returns:
            float: Oblique shock angle (radians)
        """
        maxShockAngle = ObliqueShockRelations.calc_max_shock_angle(mach, gamma)
        return brenth(
            ObliqueShockRelations.calc_theta_from_theta_beta_mach,
            maxShockAngle,
            radians(90),
            args=(mach, gamma, theta),
        )

    @staticmethod
    def calc_mach_from_theta_beta_mach(beta: float, theta: float, gamma: float) -> float:
        """ Impliments the Theta-Beta-Mach (TBM) Equation. Solves for the mach number

        Args:
            beta (float): Oblique shock angle (radians)
            theta (float): Flow deflection (wedge) angle (radians)
            gamma (float): Ratio of specific heats

        Returns:
            float: Mach number of the flow ahead of the shock wave
        """
        numerator = -2 * (1 + tan(theta) * tan(beta))
        denominator = tan(theta) * tan(beta) * (gamma + cos(2 * beta)) - 2 * (sin(beta)) ** 2
        return sqrt(numerator / denominator)

    @staticmethod
    def calc_max_flow_deflection_angle(maxShockAngle: float, mach: float, gamma: float) -> float:
        """ Calculates the max flow deflection angle for a flow

        Args:
            maxShockAngle (float): Maximum oblique shock angle (radians)
            mach (float): Mach number of flow ahead of the oblique shock
            gamma (float): Ratio of specific heats

        Returns:
            float: Mac flow deflection angle (radians)
        """
        msa = maxShockAngle
        numerator = (pow(mach, 2) * (sin(msa)) ** 2 - 1) / tan(msa)
        denominator = pow(mach, 2) * (gamma + 1) / 2 - pow(mach, 2) * (pow(sin(msa), 2)) + 1
        return atan(numerator / denominator)

    @staticmethod
    def calc_max_shock_angle(mach: float, gamma: float) -> float:
        """ Calculates the maximum oblique shock angle 

        Args:
            mach (float): Mach number of flow ahead of the shock wave
            gamma (float): Ratio of specific heats

        Returns:
            float: Maximum value of the oblique shock angle (radians)
        """
        gp1 = gamma + 1
        gm1 = gamma - 1
        # splitting up of beta_max equation
        isissq = gp1 * (1 + gm1 * pow(mach, 2) / 2 + gp1 / 16 * pow(mach, 4))
        issq = 1 / (gamma * pow(mach, 2)) * (gp1 * pow(mach, 2) / 4 + sqrt(isissq) - 1)
        return asin(sqrt(issq))

    @staticmethod
    def calc_mach_wave_angle(mach: float) -> float:
        return asin(1 / mach)

    @staticmethod
    def calc_mach_from_mach_wave_angle(machAngle: float) -> float:
        """ Calculates the Mach number fromt he mach wave angle

        Args:
            machAngle (float): Mach wave angle (mu) (radians)

        Returns:
            float: mach number
        """
        return 1 / sin(machAngle)

    def plot_theta_beta_mach_chart(self) -> None:
        """ Plots the Theta-Beta-Mach plot from the data already in the class
        """
        mach = self.mach

        machWaveAngle = degrees(ObliqueShockRelations.calc_mach_wave_angle(mach))
        maxShockAngle = degrees(ObliqueShockRelations.calc_max_shock_angle(mach, self.gamma))
        maxDeflectionAngle = degrees(
            ObliqueShockRelations.calc_max_flow_deflection_angle(radians(maxShockAngle), mach, self.gamma)
        )

        weakFlowDeflectionAngles = np.linspace(0, maxDeflectionAngle - 0.1, 100)
        weakShockAngles = np.zeros(weakFlowDeflectionAngles.shape)
        strongShockAngles = np.zeros(weakFlowDeflectionAngles.shape)

        for ii in range(weakFlowDeflectionAngles.shape[0]):
            if weakFlowDeflectionAngles[ii] == 0:
                weakShockAngles[ii] = machWaveAngle
                strongShockAngles[ii] = 90
                continue

            weakShockAngles[ii] = degrees(
                ObliqueShockRelations.calc_beta_from_theta_beta_mach_weak(
                    radians(weakFlowDeflectionAngles[ii]), mach, self.gamma
                )
            )
            strongShockAngles[ii] = degrees(
                ObliqueShockRelations.calc_beta_from_theta_beta_mach_strong(
                    radians(weakFlowDeflectionAngles[ii]), mach, self.gamma
                )
            )

        weakFlowDeflectionAngles = np.append(weakFlowDeflectionAngles, [maxDeflectionAngle])
        strongShockAngles = np.append(strongShockAngles, [maxShockAngle])
        weakShockAngles = np.append(weakShockAngles, [maxShockAngle])

        vertPointLine = np.linspace(0, self.shockAngle, 20)
        horzPointLine = np.linspace(0, self.wedgeAngle, 20)

        _, ax = plt.subplots(1, 1)

        ax.plot(weakFlowDeflectionAngles, weakShockAngles, label=" Oblique Weak Shock")
        ax.plot(weakFlowDeflectionAngles, strongShockAngles, label="Oblique Strong Shock")

        shockangle = degrees(self.shockAngle)
        wedgeAngle = degrees(self.wedgeAngle)
        ax.scatter(wedgeAngle, shockangle, label="This Flow", color="r")
        # ax.plot(np.ones(vertPointLine.shape) * self.wedgeAngle, vertPointLine, 'g')
        # ax.plot(horzPointLine, np.ones(horzPointLine.shape) * self.shockAngle, 'g')
        ax.annotate(
            "Shock {}\nWedge {}".format(round(maxShockAngle, 2), round(maxDeflectionAngle, 2)),
            xy=(maxDeflectionAngle, maxShockAngle),
            xytext=(maxDeflectionAngle + 0.02 * maxDeflectionAngle, maxShockAngle + 0.1 * maxShockAngle),
            arrowprops=dict(facecolor="black", shrink=0.05),
        )
        ax.annotate(
            "Shock {}\nWedge {}".format(round(shockangle, 2), round(wedgeAngle, 2)),
            xy=(wedgeAngle, shockangle),
            xytext=(wedgeAngle - 0.2 * wedgeAngle, shockangle - 0.4 * shockangle),
            arrowprops=dict(facecolor="black", shrink=0.05),
        )
        # ax.set_xlim(0, maxDeflectionAngle + 5)
        # ax.set_ylim(0, 90)
        ax.set_xticks([num for num in range(0, int(maxDeflectionAngle) + 5, 2)])
        ax.set_yticks([num for num in range(0, 92, 2)])
        ax.set_xlabel("Flow Deflection Angle, \u03B8 (\u00b0)")
        ax.set_ylabel("Shock Anlge, \u03B2 (\u00b0)")
        ax.set_title("Mach {} Flow".format(round(self.mach, self.precision)))
        ax.legend()
        ax.grid()
        plt.show()

