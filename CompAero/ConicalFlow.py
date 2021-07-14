from numpy import ndarray
import numpy as np
from scipy.optimize import brenth
from colorama import Back, Style, Fore
from CompAero.commonFuncs import checkValue
from CompAero.ObliqueShockRelations import ObliqueShockRelations
from GreekLetters.greekLetters import LowerCaseGreek as lcg, Misc
from CompAero.IsentropecRelations import IsentropicRelations

from typing import List, Sequence, Union
from math import isnan, sqrt, nan, pow, radians, degrees, tan, cos, sin
from dataclasses import dataclass


@dataclass
class ConicalRayRatios:
    mach: float
    tempRay_tempInit: float
    pressRay_pressInit: float
    rhoRay_rhoInit: float
    velRay_velInit: float
    cpRay: float
    rayAngle: float


class ConicalFlowRelations:
    def __init__(
        self,
        gamma: float,
        mach: float,
        increment: float = -0.01,
        coneAngle: float = nan,
        shockAngle: float = nan,
        inDegrees: bool = True,
        maxIterations: int = 10000,
    ) -> None:
        self.__gamma = gamma
        self.__mach = mach
        self.__coneAngle = coneAngle
        self.__shockAnlge = shockAngle
        self.__increment = increment
        self.__inDegrees = inDegrees
        self.__maxIter = maxIterations
        self.__velcoties = np.empty((1,))  # Store the velocites of the flow field while solving
        self.__flowAngles = np.empty((1,))  # Store the angle at which the velocites were stored

        if not self.__inDegrees:
            self.__coneAngle = degrees(self.__coneAngle)
            self.__shockAnlge = degrees(self.__shockAnlge)
        else:
            self.__increment = radians(self.__increment)

        if not checkValue(self.__gamma):
            err = "Concical flow cannot analysis cannot be completed with current value of gamma. Val: {}".format(
                self.__gamma
            )
            raise ValueError(Fore.BLACK + Back.RED + err + Style.RESET_ALL)

        if isnan(self.__increment) or self.__increment >= 0.0:
            err = "Theta increment must be a negative value. Val: {}".format(self.__increment)
            raise ValueError(Fore.BLACK + Back.RED + err + Style.RESET_ALL)

        if checkValue(self.__mach) and checkValue(self.__shockAnlge):
            self.__coneAngle = self.__calculateConeAngle(self.__shockAnlge)
            print(self.__coneAngle)

        elif checkValue(self.__mach) and checkValue(self.__coneAngle):
            self.__shockAnlge = self.__calculateShockAngle()
            print(self.__shockAnlge)

        else:
            err = "Not enough given arguments to determine conical flow field"
            raise ValueError(Fore.BLACK + Back.RED + err + Style.RESET_ALL)

        if not self.__inDegrees:
            self.__coneAngle = radians(self.__coneAngle)
            self.__shockAnlge = radians(self.__shockAnlge)
        else:
            self.__increment = degrees(self.__increment)

    @property
    def gamma(self) -> float:
        """ Ratio of specific heats """
        return self.__gamma

    @property
    def theteIncrement(self) -> float:
        """ Angle increment used when solveing Taylor-Macoll ODE"""
        return self.__increment

    @property
    def shockAngle(self) -> float:
        """ Angle of Cone shaped shock wave"""
        return self.__shockAnlge

    @property
    def coneAngle(self) -> float:
        """ Angle of the flow deflection angle (Cone) """
        return self.__coneAngle

    @property
    def mach(self) -> float:
        """ The free stream mach number seen by the cone """
        return self.__mach

    @property
    def shockJumpConditions(self) -> ObliqueShockRelations:
        return ObliqueShockRelations(
            gamma=self.__gamma, mach=self.__mach, shockAngle=self.__shockAnlge, degrees=self.__inDegrees
        )

    def __OdeFunction1(self, theta: float, Vr: float, Vtheta: float, Vr2: float, Vtheta2: float) -> float:
        return (self.__gamma - 1) / 2 * (1 - Vr2 - Vtheta2) * (2 * Vr + Vtheta / tan(theta)) - Vtheta2 * Vr

    def __OdeFunction2(self, theta: float, Vr: float, Vtheta: float) -> float:
        Vtheta2 = pow(Vtheta, 2)
        Vr2 = pow(Vr, 2)
        return (
            -1
            * self.__OdeFunction1(theta, Vr, Vtheta, Vr2, Vtheta2)
            / ((self.__gamma - 1) / 2 * (1 - Vr2 - Vtheta2) - Vtheta2)
        )

    def __initializeVelocites(self, angle: float) -> Sequence[float]:
        obliqueShockFlow = ObliqueShockRelations(self.__gamma, mach=self.__mach, shockAngle=angle)
        velocity = 1 / sqrt(2 / ((self.__gamma - 1) * pow(obliqueShockFlow.mach2, 2)) + 1)
        Vr = velocity * cos(radians(angle - obliqueShockFlow.wedgeAngle))
        Vtheta = -1 * velocity * sin(radians(angle - obliqueShockFlow.wedgeAngle))
        return [Vr, Vtheta]

    def __calculateConeAngle(
        self, shockAngle: float, offset: float = 0.0, anglesToSave: ndarray = np.array([])
    ) -> float:
        """Calculates Cone angle given shock angle and upstream mach number"""
        vr, vtheta = self.__initializeVelocites(shockAngle)
        shockAngle = radians(shockAngle)

        coneAngle = shockAngle
        angleSaveIter = 0

        currentIter = 0
        while vtheta < 0.0 and currentIter <= self.__maxIter:
            theta2 = coneAngle + self.__increment / 2
            k11 = self.__increment * vtheta
            k12 = self.__increment * self.__OdeFunction2(shockAngle, vr, vtheta)

            k21 = self.__increment * (vtheta + k12 / 2)
            k22 = self.__increment * self.__OdeFunction2(theta2, (vr + k11 / 2), (vtheta + k12 / 2))

            k31 = self.__increment * (vtheta + k22 / 2)
            k32 = self.__increment * self.__OdeFunction2(theta2, (vr + k21 / 2), (vtheta + k22 / 2))

            coneAngle += self.__increment
            k41 = self.__increment * (vtheta + k32)
            k42 = self.__increment * self.__OdeFunction2(coneAngle, (vr + k31), (vtheta + k32))

            vr += (k11 + 2 * k21 + 2 * k31 + k41) / 6
            vtheta += (k12 + 2 * k22 + 2 * k32 + k42) / 6

            currentIter += 1

            if anglesToSave.size > 0:
                if (
                    abs(np.min(np.abs(anglesToSave - round(coneAngle, 4)))) <= 1e-4
                    and angleSaveIter < self.__velcoties.shape[0]
                ):
                    self.__velcoties[angleSaveIter] = sqrt(pow(vr, 2) + pow(vtheta, 2))
                    self.__flowAngles[angleSaveIter] = coneAngle
                    angleSaveIter += 1

        return degrees(coneAngle) - offset

    def __calculateShockAngle(self) -> float:
        machWaveAngle = degrees(ObliqueShockRelations.calcMachWaveAngle(self.__mach))
        return brenth(self.__calculateConeAngle, machWaveAngle, 80, args=(self.__coneAngle,))

    def calculateConeFlowParameters(
        self, tempInit: float, gasConstR: float, angles: Union[float, ndarray, list], inDegrees: bool = True
    ) -> Union[ConicalRayRatios, dict]:
        if isinstance(angles, np.ndarray):
            self.__velcoties = np.zeros(angles.shape)
            self.__flowAngles = np.zeros(angles.shape)

        elif isinstance(angles, float) or isinstance(angles, int):
            self.__velcoties = np.array([0])
            self.__flowAngles = np.array([0])
            angles = np.array([angles])

        elif isinstance(angles, List):
            self.__velcoties = np.zeros((len(angles),))
            self.__flowAngles = np.zeros((len(angles),))
            angles = np.array(angles)

        else:
            err = "Parameter of Type {} not supported!".format(type(angles))
            raise ValueError(Fore.BLACK + Back.RED + err + Style.RESET_ALL)

        if inDegrees:
            angles = np.radians(angles)

        self.__increment = radians(self.__increment)
        self.__calculateConeAngle(self.__shockAnlge, anglesToSave=angles)

        flowParams = {}
        aInit = sqrt(self.__gamma * gasConstR * tempInit)
        shockJumpParams = ObliqueShockRelations(
            self.__gamma, mach=self.__mach, shockAngle=self.__shockAnlge, useDegrees=self.__inDegrees
        )
        conditionsInit = IsentropicRelations(self.__gamma, mach=self.__mach)

        for velocity, angle in zip(self.__velcoties, self.__flowAngles):
            params = ConicalRayRatios
            vel = sqrt(2 * self.__gamma * gasConstR / (self.__gamma - 1) * tempInit) * velocity
            mach = sqrt(
                pow(vel, 2) / (self.__gamma * gasConstR * tempInit - (self.__gamma - 1) / 2 * pow(vel, 2))
            )
            params.mach = mach
            isentropRay = IsentropicRelations(self.__gamma, mach)

            # Temperatures
            params.tempRay_tempInit = conditionsInit.t0_t / isentropRay.t0_t
            params.pressRay_pressInit = conditionsInit.p0_p * shockJumpParams.po2_po1 / isentropRay.p0_p
            params.rayAngle = degrees(angle)
            params.rhoRay_rhoInit = params.pressRay_pressInit / params.tempRay_tempInit
            params.cpRay = 2 / self.__gamma / pow(self.__mach, 2) * (params.pressRay_pressInit - 1)
            params.velRay_velInit = velocity / aInit / self.__mach
            flowParams[params.rayAngle] = params

        return flowParams

