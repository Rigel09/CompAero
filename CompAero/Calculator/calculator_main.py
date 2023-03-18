from ast import Return
from CompAero.Calculator.CalculatorUI import Ui_MainWindow
from CompAero.internal import FlowState, ShockType
from CompAero.IsentropecRelations import ISENTROPIC_VALID_OPTIONS, ISENTROPIC_CHOICE, IsentropicRelations
from CompAero.NormalShockRelations import (
    NORMAL_SHOCK_VALID_OPTIONS,
    NORMAL_SHOCK_CHOICE,
    NormalShockRelations as NSR,
)
from CompAero.ObliqueShockRelations import ObliqueShockRelations as OSR, OBLIQUE_SHOCK_VALID_OPTIONS, ObliqueShockChoice

from PyQt5.QtWidgets import QMainWindow, QWidget
from PyQt5 import QtCore, QtWidgets
import sys

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

TO_STR = lambda val, decimal: str(round(val, decimal))
PRECISION = 6


class UI(QMainWindow, Ui_MainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setupUi(self)
        self.show()
        self.layout().setContentsMargins(0, 0, 0, 0)

        # Combos
        self.isentropicOptionCombo.addItems(ISENTROPIC_VALID_OPTIONS)
        self.isentropicFlowTypeCombo.addItems([FlowState.SUPER_SONIC.name, FlowState.SUB_SONIC.name])
        self.normalShockOptionCombo.addItems(NORMAL_SHOCK_VALID_OPTIONS)
        self.obliqueShockOptionCombo.addItems(OBLIQUE_SHOCK_VALID_OPTIONS)
        self.obliqueShockTypeCombo.addItems([x.value for x in ShockType])

        # Buttons
        self.isentropicCalcBtn.clicked.connect(self.calculateIsentropicState)
        self.normalShockCalculate.clicked.connect(self.calculateNormalShockState)
        self.obliqueShockCalcBtn.clicked.connect(self.calculateObliqueShockState)
        self.obliqueShockDegreeChkBtn.clicked.connect(self.obliqueShockDegreesChkBoxUpdate)
        self.prandtlMeyerDegreeChkBtn.clicked.connect(self.prandtlMeyerDegreesChkBoxUpdate)
    
    def obliqueShockDegreesChkBoxUpdate(self) -> None:
        if self.obliqueShockDegreeChkBtn.isChecked():
            self.obliqueShockDegreeChkBtn.setText("Radians")
        else:
            self.obliqueShockDegreeChkBtn.setText("Degrees")

    def prandtlMeyerDegreesChkBoxUpdate(self) -> None:
        if self.prandtlMeyerDegreeChkBtn.isChecked():
            self.prandtlMeyerDegreeChkBtn.setText("Radians")
        else:
            self.prandtlMeyerDegreeChkBtn.setText("Degrees")

    def calculateIsentropicState(self) -> None:
        if not self.isentropicGammaEntry.text():
            return

        gamma = float(self.isentropicGammaEntry.text())

        choice = self.isentropicOptionCombo.currentText()
        if not choice:
            return

        state = None

        if choice == ISENTROPIC_CHOICE.MACH and self.isentropicMachEntry.text():
            state = IsentropicRelations(gamma=gamma, mach=float(self.isentropicMachEntry.text()))

        elif choice == ISENTROPIC_CHOICE.P0_P and self.isentropicP0PEntry.text():
            state = IsentropicRelations(gamma, p0_p=float(self.isentropicP0PEntry.text()))

        elif choice == ISENTROPIC_CHOICE.RHO0_RHO and self.isentropicRho0RhoEntry.text():
            state = IsentropicRelations(gamma, rho0_rho=float(self.isentropicRho0RhoEntry.text()))

        elif choice == ISENTROPIC_CHOICE.T0_T and self.isentropicT0TEntry.text():
            state = IsentropicRelations(gamma, t0_t=float(self.isentropicT0TEntry.text()))

        elif (
            choice == ISENTROPIC_CHOICE.A_ASTAR
            and self.isentropicAAStarEntry.text()
            and self.isentropicFlowTypeCombo.currentText()
        ):
            flowtype = FlowState(self.isentropicFlowTypeCombo.currentText())
            aaStar = float(self.isentropicAAStarEntry.text())
            state = IsentropicRelations(gamma, a_aStar=aaStar, flowType=flowtype)

        if state is not None:
            self.isentropicMachEntry.setText(TO_STR(state.mach, PRECISION))
            self.isentropicP0PEntry.setText(TO_STR(state.p0_p, PRECISION))
            self.isentropicAAStarEntry.setText(TO_STR(state.a_aStar, PRECISION))
            self.isentropicT0TEntry.setText(TO_STR(state.t0_t, PRECISION))
            self.isentropicRho0RhoEntry.setText(TO_STR(state.rho0_rho, PRECISION))

    def calculateNormalShockState(self) -> None:
        if not self.normalShockGammaEntry.text():
            return

        gamma = float(self.normalShockGammaEntry.text())

        choice = NORMAL_SHOCK_CHOICE(self.normalShockOptionCombo.currentText())
        if not choice:
            return
        state = None

        if choice == NORMAL_SHOCK_CHOICE.MACH and self.normalShockM1Entry.text():
            state = NSR(gamma, mach=float(self.normalShockM1Entry.text()))

        elif choice == NORMAL_SHOCK_CHOICE.P2_P1 and self.normalShockP2P1Entry.text():
            state = NSR(gamma, p2_p1=float(self.normalShockP2P1Entry.text()))

        elif choice == NORMAL_SHOCK_CHOICE.RHO2_RHO1 and self.normalShockRho2Rho1Entry.text():
            state = NSR(gamma, rho2_rho1=float(self.normalShockRho2Rho1Entry.text()))

        elif choice == NORMAL_SHOCK_CHOICE.T2_T1 and self.normalShockT2T1Entry.text():
            state = NSR(gamma, t2_t1=float(self.normalShockT2T1Entry.text()))

        elif choice == NORMAL_SHOCK_CHOICE.PO2_PO1 and self.normalShockP02P01Entry.text():
            state = NSR(gamma, po2_po1=float(self.normalShockP02P01Entry.text()))

        elif choice == NORMAL_SHOCK_CHOICE.PO2_P1 and self.normalShockP02P1Entry.text():
            state = NSR(gamma, po2_p1=float(self.normalShockP02P1Entry.text()))

        elif choice == NORMAL_SHOCK_CHOICE.M2 and self.normalShockM2Entry.text():
            state = NSR(gamma, m2=float(self.normalShockM2Entry.text()))

        else:
            print("Invalid Choice ")

        if state is not None:
            self.normalShockM1Entry.setText(TO_STR(state.mach, PRECISION))
            self.normalShockM2Entry.setText(TO_STR(state.mach2, PRECISION))
            self.normalShockP2P1Entry.setText(TO_STR(state.p2_p1, PRECISION))
            self.normalShockRho2Rho1Entry.setText(TO_STR(state.rho2_rho1, PRECISION))
            self.normalShockT2T1Entry.setText(TO_STR(state.t2_t1, PRECISION))
            self.normalShockP02P01Entry.setText(TO_STR(state.po2_po1, PRECISION))
            self.normalShockP02P1Entry.setText(TO_STR(state.po2_p1, PRECISION))

    def calculateObliqueShockState(self) -> None:
        if not self.obliqueShockGammaEntry.text():
            return
        
        gamma = float(self.obliqueShockGammaEntry.text())
        choice = ObliqueShockChoice(self.obliqueShockOptionCombo.currentText())
        shockType = ShockType(self.obliqueShockTypeCombo.currentText())
        wedgeAngleValid = bool(self.obliqueShockWedgeAngleEdit.text())
        shockAngleValid = bool(self.obliqueShockAngleEdit.text())
        m1Valid         = bool(self.obliqueShockM1Edit.text())
        mn1Valid        = bool(self.obliqueShockMn1Edit.text())
        mn2Valid        = bool(self.obiqueShockMn2Edit.text())
        m2Valid         = bool(self.obliqueShockM2Edit.text())
        p2p1Valid       = bool(self.obliqueShockP2P1Edit.text())
        rho2rho1Valid   = bool(self.obliqueShockRho2Rho1Edit.text())
        po2p1Valid      = bool(self.obliqueShockPo2P1Edit.text())
        t2t1Valid       = bool(self.obliqueShockT2T1Edit.text())
        po2po1Valid     = bool(self.obliqueShockPo2Po1Edit.text())
        degrees         = bool(not self.obliqueShockDegreeChkBtn.isChecked())
        
        if not choice:
            return
        
        state = None
        
        # common options
        opts = {
            "gamma": gamma,
            "useDegrees": degrees,
            "shockType": shockType
        }
        
        if choice == ObliqueShockChoice.MACH_WEDGE_ANGLE and m1Valid and wedgeAngleValid:
            m = float(self.obliqueShockM1Edit.text())
            wa = float(self.obliqueShockWedgeAngleEdit.text())
            state = OSR(**opts, mach=m, wedgeAngle=wa)
            
        if not shockAngleValid:
            return
        
        opts["shockAngle"] = float(self.obliqueShockAngleEdit.text())

        if choice == ObliqueShockChoice.MACH_SHOCK_ANGLE and m1Valid:
            m = float(self.obliqueShockM1Edit.text())
            state = OSR(**opts, mach=m)
        
        elif choice == ObliqueShockChoice.MACH_N_1_SHOCK_ANGLE and mn1Valid:
            m = float(self.obliqueShockMn1Edit.text())
            state = OSR(**opts, mn1=m)
        
        elif choice == ObliqueShockChoice.M2_WEDGE_SHOCK_ANGLE and m2Valid and wedgeAngleValid:
            m2 = float(self.obliqueShockM2Edit.text())
            wa = float(self.obliqueShockWedgeAngleEdit.text())
            state = OSR(**opts, m2=m2, wedgeAngle=wa)
        
        elif choice == ObliqueShockChoice.P2_P1_SHOCK_ANGLE and p2p1Valid:
            p2p1 = float(self.obliqueShockP2P1Edit.text())
            state = OSR(**opts, p2_p1=p2p1)
        
        elif choice == ObliqueShockChoice.RHO2_RHO1_SHOCK_ANGLE and rho2rho1Valid:
            r2r1 = float(self.obliqueShockRho2Rho1Edit.text())
            state = OSR(**opts, rho2_rho1=r2r1)
        
        elif choice == ObliqueShockChoice.T2_T1_SHOCK_ANGLE and t2t1Valid:
            t2t1 = float(self.obliqueShockT2T1Edit.text())
            state = OSR(**opts, t2_t1=t2t1)
        
        elif choice == ObliqueShockChoice.PO2_PO1_SHOCK_ANGLE and po2po1Valid:
            p = float(self.obliqueShockPo2Po1Edit.text())
            state = OSR(**opts, po2_po1=p)
        
        elif choice == ObliqueShockChoice.PO2_P1_SHOCK_ANGLE and po2p1Valid:
            p = float(self.obliqueShockPo2P1Edit.text())
            state = OSR(**opts, po2_p1=p)
        
        elif choice == ObliqueShockChoice.MN2_SHOCK_ANGLE and mn2Valid:
            m = float(self.obiqueShockMn2Edit.text())
            state = OSR(**opts, mn2=m)
        
        if state:
            self.obliqueShockWedgeAngleEdit.setText(TO_STR(state.wedgeAngle, PRECISION))
            self.obliqueShockAngleEdit.setText(TO_STR(state.shockAngle, PRECISION))
            self.obliqueShockM1Edit.setText(TO_STR(state.mach, PRECISION))
            self.obliqueShockMn1Edit.setText(TO_STR(state.machNorm1, PRECISION))
            self.obiqueShockMn2Edit.setText(TO_STR(state.machNorm2, PRECISION))
            self.obliqueShockM2Edit.setText(TO_STR(state.mach2, PRECISION))
            self.obliqueShockP2P1Edit.setText(TO_STR(state.p2_p1, PRECISION))
            self.obliqueShockRho2Rho1Edit.setText(TO_STR(state.rho2_rho1, PRECISION))
            self.obliqueShockPo2P1Edit.setText(TO_STR(state.po2_p1, PRECISION))
            self.obliqueShockT2T1Edit.setText(TO_STR(state.t2_t1, PRECISION))
            self.obliqueShockPo2Po1Edit.setText(TO_STR(state.po2_po1, PRECISION))

    def calculateFannoFlowState(self) -> None:
        pass

    def calculateRayleighFlowState(self) -> None:
        pass

    def calculatePrandtlMeyerState(self) -> None:
        pass

    def calculateRocketNozzleState(self) -> None:
        pass


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    _ = UI()
    app.exec_()


if __name__ == "__main__":
    main()