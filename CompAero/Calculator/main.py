from ast import Return
from CompAero.Calculator.CalculatorUI import Ui_MainWindow
from CompAero.internal import FlowState, ShockType
from CompAero.IsentropecRelations import ISENTROPIC_VALID_OPTIONS, ISENTROPIC_CHOICE, IsentropicRelations

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

        self.isentropicOptionCombo.addItems(ISENTROPIC_VALID_OPTIONS)
        self.isentropicCalcBtn.clicked.connect(self.calculateIsentropicState)
        self.isentropicFlowTypeCombo.addItems([FlowState.SUPER_SONIC.name, FlowState.SUB_SONIC.name])

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


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    _ = UI()
    app.exec_()
