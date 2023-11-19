"""
This file contains the relationships for oblique shock waves
"""

from enum import Enum
from math import asin, atan, cos, degrees, nan, radians, sin, sqrt, tan

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from scipy.optimize import brenth  # type: ignore

from CompAero.greek_letters import LowerCaseGreek as lcg
from CompAero.internal import (
    GammaNotDefinedError,
    InvalidOptionCombinationError,
    ShockType,
    check_value,
    footer,
    named_header,
    named_subheader,
    seperator,
    to_string,
)
from CompAero.normal_shock_relations import NormalShockRelations


class ObliqueShockChoice(Enum):
    """Valid Oblique Shock Initialization Choices"""

    MACH_SHOCK_ANGLE = "gamma, mach, shock angle"
    MACH_WEDGE_ANGLE = "gamma, mach, wedge angle"
    MACH_N_1_SHOCK_ANGLE = "gamma, mach normal 1, shock angle "
    M2_WEDGE_SHOCK_ANGLE = "gamma, mach behind shock (M2), wedge angle, shock angle"
    P2_P1_SHOCK_ANGLE = "gamma, P2/P1, shock angle"
    RHO2_RHO1_SHOCK_ANGLE = "gamma, Rho2/Rho1, shock angle"
    T2_T1_SHOCK_ANGLE = "gamma, T2/T1, shock angle"
    PO2_PO1_SHOCK_ANGLE = "gamma, P02/P01, shock angle"
    PO2_P1_SHOCK_ANGLE = "gamma, P02/P1, shock angle"
    MN2_SHOCK_ANGLE = "gamma, MN2, shock angle"


OBLIQUE_SHOCK_VALID_OPTIONS = [x.value for x in ObliqueShockChoice]


class ObliqueShockRelations(NormalShockRelations):
    """
    This class is a collective name space for basic calculations regarding Oblique Shock Properties
    The constructor of this class can  determine the entire state of the flow given a partial state

    Args:
        gamma (float): Ratio of specific heats
        use_degrees (bool, optional): If true then the angles provided to the class are in degrees.
                                        Defaults to True.
        shock_angle (float, optional): Angle of the oblique shock wave. Defaults to nan.
        wedge_angle (float, optional): Angle of the wedge that is deflecting the flow. Defaults to
                                        nan.
        mn1 (float, optional): Normal component of mach number ahead of the shock. Defaults to nan.
        mn2 (float, optional): Normal component of mach number behind the shock. Defaults to nan.
        mach (float, optional): Mach number ahead of the shock. Defaults to nan.
        p2_p1 (float, optional): Ratio of pressure behind shock to ahead of shock. Defaults to nan.
        rho2_rho1 (float, optional): Ratio of density behind shock to ahead of shock. Defaults to
                                        nan.
        t2_t1 (float, optional): Ratio of temp behind shock to ahead of shock. Defaults to nan.
        po2_po1 (float, optional): Ratio of total pressure to behind shock to ahead of shock.
                                    Defaults to nan.
        po2_p1 (float, optional): Ratio Total pressure behind shock to static pressure ahead of
                                    shock. Defaults to nan.
        m2 (float, optional): mach number behind shock wave. Defaults to nan.
        shock_type (ShockType, optional): Describes whether the oblique shock is Weak or Strong.
                                            Defaults to ShockType.WEAK.

    Raises:
        RuntimeError: If shock wave angle is found to be less than the mach wave angle at
        RuntimeError: If the deflection angle is greater than the max deflection angle for the flow
        GammaNotDefinedError: Raised if Gamma is undefined
        InvalidOptionCombinationError: Raised if an invalid combination of parameters is given and
                                        the flow state cannot be determined

    Valid_Combinations_of_Parameters:
        1: gamma, mach, shock angle\n
        2. gamma, mach, wedge angle\n
        2: gamma, mach normal 1, shock angle \n
        3: gamma, mach behind shock, wedge angle, shock angle\n
        4: gamma, shock angle, wedge angle\n
        5: gamma, P2/P1, shock angle\n
        6: gamma, Rho2/Rho1, shock angle\n
        7: gamma, T2/T1, shock angle\n
        8: gamma, P02/P01, shock angle\n
        9: gamma, P02/P1, shock angle\n
        10: gamma, dwn_strm_mach, shock angle\n
    """

    # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        gamma: float,
        use_degrees: bool = True,
        shock_angle: float = nan,
        wedge_angle: float = nan,
        mn1: float = nan,
        mn2: float = nan,
        mach: float = nan,
        p2_p1: float = nan,
        rho2_rho1: float = nan,
        t2_t1: float = nan,
        po2_po1: float = nan,
        po2_p1: float = nan,
        m2: float = nan,
        shock_type: ShockType = ShockType.WEAK,
    ) -> None:
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-statements
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-locals
        self.use_degrees = use_degrees
        """ Angles passed to the constructor are in degrees"""
        self.shock_angle = shock_angle
        """ Angle Oblique shock"""
        self.wedge_angle = wedge_angle
        """ Angle of the flow deflection (wedge) angle"""
        self.mach_normal_1 = mn1
        """ The normal component of the mach number ahead the shock wave"""
        self.mach_normal_2 = mn2
        """ The normal component of the mach number behind the shock wave"""
        self.gamma = gamma
        """ Ratio of specific heats """
        self.mach = mach
        """ Mach number of the flow """
        self.p2_p1 = p2_p1
        """ Ratio of pressure behind the shock wave to pressure before the shock wave P2/P1 """
        self.rho2_rho1 = rho2_rho1
        """ Ratio of density behind the shock to density before the shock rho2/rho1 """
        self.t2_t1 = t2_t1
        """ Ratio of temperature behind the shock to temperature before the shock wave T2/T1 """
        self.po2_po1 = po2_po1
        """ Ratio of pressure behind the shock wave to pressure before the shock wave P2/P1 """
        self.po2_p1 = po2_p1
        """ Ratio of total pressure behind the shock to pressure before the shock P02/P1 """
        self.mach2 = m2
        """ Mach number behind the shock wave """
        self.shock_type: ShockType = shock_type
        """ The type of shock whether it is strong or weak """
        self.precision = 4
        """ Precision to use when printing output to the console defaults to four """

        if self.use_degrees and check_value(self.shock_angle):
            self.shock_angle = radians(self.shock_angle)
        if self.use_degrees and check_value(self.wedge_angle):
            self.wedge_angle = radians(self.wedge_angle)

        if not check_value(self.gamma):
            raise GammaNotDefinedError()

        if check_value(self.mach, self.shock_angle):
            mach_wave_angle = ObliqueShockRelations.calc_mach_wave_angle(self.mach)
            if self.shock_angle < mach_wave_angle or self.shock_angle > radians(90):
                sa = degrees(self.shock_angle)
                mwa = degrees(mach_wave_angle)
                m = self.mach
                raise ValueError(
                    f"Shock Angle of [{sa}] must be between mach wave angle of [{mwa}] and 90"
                    f" degrees for mach [{m}]"
                )

        elif check_value(self.mach, self.wedge_angle):
            max_shock_angle = ObliqueShockRelations.calc_max_shock_angle(self.mach, self.gamma)
            max_wedge_angle = ObliqueShockRelations.calc_max_flow_deflection_angle(
                max_shock_angle, self.mach, self.gamma
            )
            if self.wedge_angle > max_wedge_angle:
                wa = degrees(self.wedge_angle)
                mwa = degrees(max_wedge_angle)
                m = self.mach
                raise ValueError(
                    f"Wedge angle of [{wa}] is greater than a maximum wedge angle of [{mwa}] for"
                    f"mach [{m}]"
                )

        elif check_value(self.mach_normal_1, self.shock_angle):
            if self.mach_normal_1 < 1:
                raise ValueError(
                    "Normal Component of mach ahead of shock wave must be greater than 1"
                )
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(m2, self.wedge_angle, self.shock_angle):
            self.mach_normal_2 = (
                ObliqueShockRelations.calc_mach_normal_behind_shock_from_mach_behind_shock(
                    m2, self.shock_angle, self.wedge_angle
                )
            )
            if self.mach_normal_2 > 1:
                raise ValueError("Normal component of downstream mach number has to be less than 1")
            self.mach_normal_1 = (
                NormalShockRelations.calc_mach_before_normal_shock_from_mach_after_shock(
                    self.mach_normal_2, self.gamma
                )
            )
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(self.shock_angle, self.wedge_angle):
            self.mach = ObliqueShockRelations.calc_mach_from_theta_beta_mach(
                self.shock_angle, self.wedge_angle, self.gamma
            )
            mach_angle = ObliqueShockRelations.calc_mach_wave_angle(self.mach)
            max_shock_angle = ObliqueShockRelations.calc_max_shock_angle(self.mach, self.gamma)
            max_deflection_angle = ObliqueShockRelations.calc_max_flow_deflection_angle(
                max_shock_angle, self.mach, self.gamma
            )

            if mach_angle > self.shock_angle:
                sa = round(self.shock_angle, 4)
                ma = round(mach_angle, 4)
                raise RuntimeError(
                    f"Shock Angle, {sa}, cannot be less than that of Mach Wave, {ma}"
                )

            if max_deflection_angle < self.wedge_angle:
                wa = round(self.wedge_angle, 4)
                mda = round(max_deflection_angle, 4)
                raise RuntimeError(f"SHOCK DETACHED!! Deflection Angles: Given: {wa}\tMax: {mda}")

        elif check_value(po2_p1, self.shock_angle):
            self.mach_normal_1 = NormalShockRelations.calc_mach_from_po2_p1(po2_p1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(p2_p1, self.shock_angle):
            self.mach_normal_1 = NormalShockRelations.calc_mach_from_p2_p1(p2_p1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(rho2_rho1, self.shock_angle):
            self.mach_normal_1 = NormalShockRelations.calc_mach_from_rho2_rho1(
                rho2_rho1, self.gamma
            )
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(t2_t1, self.shock_angle):
            self.mach_normal_1 = NormalShockRelations.calc_mach_from_t2_t1(t2_t1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(po2_po1, self.shock_angle):
            self.mach_normal_1 = NormalShockRelations.calc_mach_from_po2_po1(po2_po1, self.gamma)
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        elif check_value(self.mach_normal_2, self.shock_angle):
            self.mach_normal_1 = (
                NormalShockRelations.calc_mach_before_normal_shock_from_mach_after_shock(
                    self.mach_normal_2, self.gamma
                )
            )
            self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
                self.mach_normal_1, self.shock_angle
            )

        else:
            raise InvalidOptionCombinationError()

        if check_value(self.mach) and (
            (check_value(self.shock_angle) or check_value(self.wedge_angle))
        ):
            self._calc_state()

        else:
            raise RuntimeError(
                "Values for mach and shock angle or wedge angle are needed to"
                "determine flow state. For some reason those wasn't computed."
            )

        super().__init__(self.gamma, mach=self.mach_normal_1)
        self.mach = ObliqueShockRelations.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
            self.mach_normal_1, self.shock_angle
        )
        self.mach_normal_2 = self.mach2
        self.mach2 = ObliqueShockRelations.calc_mach_behind_shock(
            self.mach_normal_2, self.wedge_angle, self.shock_angle
        )

        print(f"Here use degrees: {self.use_degrees} {self.wedge_angle}  {self.shock_angle}")

        if self.use_degrees:
            self.wedge_angle = degrees(self.wedge_angle)
            self.shock_angle = degrees(self.shock_angle)
            print(f"Here use degrees: {self.use_degrees} {self.wedge_angle}  {self.shock_angle}")

    def _calc_state(self) -> None:
        if check_value(self.wedge_angle):
            if self.shock_type == ShockType.WEAK:
                self.shock_angle = ObliqueShockRelations.calc_beta_from_theta_beta_mach_weak(
                    self.wedge_angle, self.mach, self.gamma
                )

            elif self.shock_type == ShockType.STRONG:
                self.shock_angle = ObliqueShockRelations.calc_beta_from_theta_beta_mach_strong(
                    self.wedge_angle, self.mach, self.gamma
                )

            else:
                raise ValueError(
                    f"Incorrect shock type specified -> [{self.shock_type}]. Choices: Strong, Weak"
                )

        elif check_value(self.shock_angle):
            self.wedge_angle = ObliqueShockRelations.calc_theta_from_theta_beta_mach(
                self.shock_angle, self.mach, self.gamma
            )

        self.mach_normal_1 = ObliqueShockRelations.calc_mach_normal_ahead_shock(
            self.mach, self.shock_angle
        )

    def __str__(self) -> str:
        return "".join(
            [
                named_header("Oblique Shock Relations at Mach", self.mach, self.precision),
                seperator(),
                named_subheader("Upstream Conditions"),
                to_string(lcg.gamma, self.gamma, self.precision),
                to_string("Mach", self.mach, self.precision, dot_line=True),
                to_string("Mach Normal Component", self.mach_normal_1, self.precision),
                to_string(
                    f"Flow Deflection Angle {lcg.theta}",
                    self.wedge_angle,
                    self.precision,
                    dot_line=True,
                ),
                to_string(f"Shock Angle {lcg.beta}", self.shock_angle, self.precision),
                to_string("Flow Turn Type", self.shock_type.name, self.precision),
                seperator(),
                named_subheader("Shock Jump Conditions"),
                to_string("P2/P1", self.p2_p1, self.precision),
                to_string(
                    "{}2/{}1".format(*[lcg.rho] * 2), self.rho2_rho1, self.precision, dot_line=True
                ),
                to_string("T2/T1", self.t2_t1, self.precision),
                to_string("P02/P01", self.po2_po1, self.precision, dot_line=True),
                to_string("P02/P1", self.po2_p1, self.precision),
                seperator(),
                named_subheader("Downstream Conditions"),
                to_string("Mach", self.mach2, self.precision, dot_line=True),
                to_string("Mach Normal Component", self.mach_normal_2, self.precision),
                footer(),
            ]
        )

    @staticmethod
    def calc_mach_normal_ahead_shock(mach: float, beta: float) -> float:
        """Calculates the normal component of the mach number ahead of the shock wave

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

        mach_wave = ObliqueShockRelations.calc_mach_wave_angle(mach)

        if abs(beta - mach_wave) < 1e-5:
            return 1.0

        if beta < mach_wave:
            return nan

        return mach * sin(beta)

    @staticmethod
    def calc_mach_ahead_shock_from_mach_normal_ahead_shock(
        mach_normal_1: float, beta: float
    ) -> float:
        """Calculates the upstream mach number from the normal component of the upstream mach number

        Args:
            machNormal1 (float): Normal Component of mach number ahead of the shock wave
            beta (float): Angle of oblique shock in radians

        Returns:
            float: Returns value of mach number of the flow ahead of the shock wave
        """
        return mach_normal_1 / sin(beta)

    @staticmethod
    def calc_beta_from_mach_mach_normal_ahead_shock(mach: float, mach_normal_1: float) -> float:
        """
        Calculates the Oblique shock angle from the normal component of the mach number that is
        ahead of the shock wave

        Args:
            mach (float): Mach number of the flow ahead of the shock wave
            mach_normal_1  (float): Normal Component of mach number of the flow that is ahead of an
                                    oblique shock wave

        Returns:
            float: Oblique Shock Angle
        """
        return asin(mach_normal_1 / mach)

    @staticmethod
    def calc_mach_behind_shock(mach_normal_2: float, theta: float, beta: float) -> float:
        """Calculates the Mach number behind the oblique shock wave

        Args:
            mach_normal_2 (float): Normal Component of the mach number behind the shock wave
            theta (float): Flow Deflection Angle (radians) (Wedge angle)
            beta (float): Oblique shock angle (radians)

        Returns:
            float: Mach number of flow behind the oblique shock
        """
        return mach_normal_2 / sin(beta - theta)

    @staticmethod
    def calc_mach_normal_behind_shock_from_mach_behind_shock(
        mach2: float, beta: float, theta: float
    ) -> float:
        """Calculates the normal component of the mach number behind the oblique shock

        Args:
            mach2 (float): Mach number of flow behind the oblique shock
            beta (float): Oblique shock angle (radians)
            theta (float): Flow deflections (Wedge) angle (radians)

        Returns:
            float: Normal component of mach numnber of the flow behind the shock wave
        """
        return mach2 * sin(beta - theta)

    @staticmethod
    def calc_theta_from_theta_beta_mach(
        beta: float, mach: float, gamma: float, offset: float = 0.0
    ) -> float:
        """Impliments the Theta-Beta-Mach (TBM) equation. Solves for Theta

        Args:
            beta (float): Oblique shock angle (radians)
            mach (float): Mach number of flow ahead of the shock wave
            gamma (float): Ratio of specific heats
            offset (float, optional): [description]. Defaults to 0.0.

        Returns:
            float: Flow deflection (Wedge) angle (radians)
        """
        m_sqr = pow(mach, 2)
        num = m_sqr * pow(sin(beta), 2) - 1
        denom = m_sqr * (gamma + cos(2 * beta)) + 2
        theta = atan(2 * 1 / tan(beta) * num / denom) - offset
        return theta

    @staticmethod
    def calc_beta_from_theta_beta_mach_weak(theta: float, mach: float, gamma: float) -> float:
        """
        Impliments the Theta-Beta-Mach (TBM) equation. Solves for Beta (shock angle) assuming the
        shock is weak

        Args:
            theta (float): Flow deflection (Wedge) angle (radians)
            mach (float): Mach number of the flow ahead of the oblique shock
            gamma (float): ratio of specific heats

        Returns:
            float: Oblique shock angle (radians)
        """
        max_shock_angle = ObliqueShockRelations.calc_max_shock_angle(mach, gamma)
        min_shock_angle = ObliqueShockRelations.calc_mach_wave_angle(mach)
        return brenth(
            ObliqueShockRelations.calc_theta_from_theta_beta_mach,
            min_shock_angle,
            max_shock_angle,
            args=(mach, gamma, theta),
        )  # type: ignore

    @staticmethod
    def calc_beta_from_theta_beta_mach_strong(theta: float, mach: float, gamma: float) -> float:
        """
        Impliments the Theta-Beta-Mach (TBM) equation. Solves for Beta (shock angle) assuming a
        strong shock wave

        Args:
            theta (float): Flow deflection (Wedge) angle (radians)
            mach (float): Mach number of the flow ahead of the oblique shock
            gamma (float): ratio of specific heats

        Returns:
            float: Oblique shock angle (radians)
        """
        maxshock_angle = ObliqueShockRelations.calc_max_shock_angle(mach, gamma)
        return brenth(
            ObliqueShockRelations.calc_theta_from_theta_beta_mach,
            maxshock_angle,
            radians(90),
            args=(mach, gamma, theta),
        )  # type: ignore

    @staticmethod
    def calc_mach_from_theta_beta_mach(beta: float, theta: float, gamma: float) -> float:
        """Impliments the Theta-Beta-Mach (TBM) Equation. Solves for the mach number

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
    def calc_max_flow_deflection_angle(maxshock_angle: float, mach: float, gamma: float) -> float:
        """Calculates the max flow deflection angle for a flow

        Args:
            maxshock_angle (float): Maximum oblique shock angle (radians)
            mach (float): Mach number of flow ahead of the oblique shock
            gamma (float): Ratio of specific heats

        Returns:
            float: Mac flow deflection angle (radians)
        """
        msa = maxshock_angle
        numerator = (pow(mach, 2) * (sin(msa)) ** 2 - 1) / tan(msa)
        denominator = pow(mach, 2) * (gamma + 1) / 2 - pow(mach, 2) * (pow(sin(msa), 2)) + 1
        return atan(numerator / denominator)

    @staticmethod
    def calc_max_shock_angle(mach: float, gamma: float) -> float:
        """Calculates the maximum oblique shock angle

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
        """
        Calculates the mach wave angle given the mach number

        Args:
            mach float: Mach number

        Returns (float): Returns the mach wave angle in radians
        """
        return asin(1 / mach)

    @staticmethod
    def calc_mach_from_mach_wave_angle(mach_angle: float) -> float:
        """Calculates the Mach number fromt he mach wave angle

        Args:
            mach_angle (float): Mach wave angle (mu) (radians)

        Returns:
            float: mach number
        """
        return 1 / sin(mach_angle)

    def plot_theta_beta_mach_chart(self) -> None:
        """Plots the Theta-Beta-Mach plot from the data already in the class"""
        mach = self.mach

        mach_wave_angle = degrees(ObliqueShockRelations.calc_mach_wave_angle(mach))
        max_shock_angle = degrees(ObliqueShockRelations.calc_max_shock_angle(mach, self.gamma))
        max_deflection_angle = degrees(
            ObliqueShockRelations.calc_max_flow_deflection_angle(
                radians(max_shock_angle), mach, self.gamma
            )
        )

        weak_deflection_angles = np.linspace(0, max_deflection_angle - 0.1, 100)
        weak_shock_angles = np.zeros(weak_deflection_angles.shape)
        strong_shock_angles = np.zeros(weak_deflection_angles.shape)

        for ii in range(weak_deflection_angles.shape[0]):
            if weak_deflection_angles[ii] == 0:
                weak_shock_angles[ii] = mach_wave_angle
                strong_shock_angles[ii] = 90
                continue

            weak_shock_angles[ii] = degrees(
                ObliqueShockRelations.calc_beta_from_theta_beta_mach_weak(
                    radians(weak_deflection_angles[ii]), mach, self.gamma
                )
            )
            strong_shock_angles[ii] = degrees(
                ObliqueShockRelations.calc_beta_from_theta_beta_mach_strong(
                    radians(weak_deflection_angles[ii]), mach, self.gamma
                )
            )

        weak_deflection_angles = np.append(weak_deflection_angles, [max_deflection_angle])
        strong_shock_angles = np.append(strong_shock_angles, [max_shock_angle])
        weak_shock_angles = np.append(weak_shock_angles, [max_shock_angle])

        # vertPointLine = np.linspace(0, self.shock_angle, 20)
        # horzPointLine = np.linspace(0, self.wedge_angle, 20)

        _, ax = plt.subplots(1, 1)

        ax.plot(weak_deflection_angles, weak_shock_angles, label=" Oblique Weak Shock")
        ax.plot(weak_deflection_angles, strong_shock_angles, label="Oblique Strong Shock")

        shock_angle = degrees(self.shock_angle)
        wedge_angle = degrees(self.wedge_angle)
        ax.scatter(wedge_angle, shock_angle, label="This Flow", color="r")
        # ax.plot(np.ones(vertPointLine.shape) * self.wedge_angle, vertPointLine, 'g')
        # ax.plot(horzPointLine, np.ones(horzPointLine.shape) * self.shock_angle, 'g')
        ax.annotate(
            f"Shock {round(max_shock_angle, 2)}\nWedge {round(max_deflection_angle, 2)}",
            xy=(max_deflection_angle, max_shock_angle),
            xytext=(
                max_deflection_angle + 0.02 * max_deflection_angle,
                max_shock_angle + 0.1 * max_shock_angle,
            ),
            arrowprops={"facecolor": "black", "shrink": 0.05},
        )
        ax.annotate(
            f"Shock {round(shock_angle, 2)}\nWedge {round(wedge_angle, 2)}",
            xy=(wedge_angle, shock_angle),
            xytext=(wedge_angle - 0.2 * wedge_angle, shock_angle - 0.4 * shock_angle),
            arrowprops={"facecolor": "black", "shrink": 0.05},
        )
        # ax.set_xlim(0, maxDeflectionAngle + 5)
        # ax.set_ylim(0, 90)
        ax.set_xticks(list(range(0, int(max_deflection_angle) + 5, 2)))
        ax.set_yticks(list(range(0, 92, 2)))
        ax.set_xlabel("Flow Deflection Angle, \u03B8 (\u00b0)")
        ax.set_ylabel("Shock Anlge, \u03B2 (\u00b0)")
        ax.set_title(f"Mach {round(self.mach, self.precision)} Flow")
        ax.legend()
        ax.grid()
        plt.show()
