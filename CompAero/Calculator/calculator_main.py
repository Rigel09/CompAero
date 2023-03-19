from typing import Optional
from CompAero.Calculator.CalculatorUI import Ui_MainWindow
from CompAero.internal import FlowState, ShockType
from CompAero.IsentropecRelations import ISENTROPIC_VALID_OPTIONS, ISENTROPIC_CHOICE, IsentropicRelations
from CompAero.NormalShockRelations import (
    NORMAL_SHOCK_VALID_OPTIONS,
    NORMAL_SHOCK_CHOICE,
    NormalShockRelations as NSR,
)
from CompAero.ObliqueShockRelations import (
    ObliqueShockRelations as OSR,
    OBLIQUE_SHOCK_VALID_OPTIONS,
    ObliqueShockChoice
)
from CompAero.FannoFlowRelations import (
    FannoFlowRelations as FFR,
    FANNO_FLOW_VALID_OPTIONS,
    FannoFlowChoice
)
from PyQt5.QtWidgets import QMainWindow, QWidget
from PyQt5 import QtCore, QtWidgets
import sys

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

PRECISION = 6
TO_STR = lambda val, decimal = PRECISION: str(round(val, decimal))


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
        self.fannoOptionCombo.addItems(FANNO_FLOW_VALID_OPTIONS)
        self.fannoFlowTypeCombo.addItems([x.value for x in FlowState])

        # Buttons
        self.isentropicCalcBtn.clicked.connect(self.calculateIsentropicState)
        self.normalShockCalculate.clicked.connect(self.calculateNormalShockState)
        self.obliqueShockCalcBtn.clicked.connect(self.calculateObliqueShockState)
        self.obliqueShockDegreeChkBtn.clicked.connect(self.obliqueShockDegreesChkBoxUpdate)
        self.prandtlMeyerDegreeChkBtn.clicked.connect(self.prandtlMeyerDegreesChkBoxUpdate)
        self.fannoCalculateBtn.clicked.connect(self.calculateFannoFlowState)
        self.fannoApplyPipeParamBtn.clicked.connect(self.calculateFannoFrictionAddition)
        
        # Saved states
        self._fannoState: Optional[FFR] = None
    
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
            self.isentropicMachEntry.setText(TO_STR(state.mach))
            self.isentropicP0PEntry.setText(TO_STR(state.p0_p))
            self.isentropicAAStarEntry.setText(TO_STR(state.a_aStar))
            self.isentropicT0TEntry.setText(TO_STR(state.t0_t))
            self.isentropicRho0RhoEntry.setText(TO_STR(state.rho0_rho))

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
            self.normalShockM1Entry.setText(TO_STR(state.mach))
            self.normalShockM2Entry.setText(TO_STR(state.mach2))
            self.normalShockP2P1Entry.setText(TO_STR(state.p2_p1))
            self.normalShockRho2Rho1Entry.setText(TO_STR(state.rho2_rho1))
            self.normalShockT2T1Entry.setText(TO_STR(state.t2_t1))
            self.normalShockP02P01Entry.setText(TO_STR(state.po2_po1))
            self.normalShockP02P1Entry.setText(TO_STR(state.po2_p1))

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
            self.obliqueShockWedgeAngleEdit.setText(TO_STR(state.wedgeAngle))
            self.obliqueShockAngleEdit.setText(TO_STR(state.shockAngle))
            self.obliqueShockM1Edit.setText(TO_STR(state.mach))
            self.obliqueShockMn1Edit.setText(TO_STR(state.machNorm1))
            self.obiqueShockMn2Edit.setText(TO_STR(state.machNorm2))
            self.obliqueShockM2Edit.setText(TO_STR(state.mach2))
            self.obliqueShockP2P1Edit.setText(TO_STR(state.p2_p1))
            self.obliqueShockRho2Rho1Edit.setText(TO_STR(state.rho2_rho1))
            self.obliqueShockPo2P1Edit.setText(TO_STR(state.po2_p1))
            self.obliqueShockT2T1Edit.setText(TO_STR(state.t2_t1))
            self.obliqueShockPo2Po1Edit.setText(TO_STR(state.po2_po1))

    def calculateFannoFlowState(self) -> None:
        if not self.fannoGammaEntry.text():
            return
        
        gamma = float(self.fannoGammaEntry.text())
        choice = FannoFlowChoice(self.fannoOptionCombo.currentText())
        shockType = FlowState(self.fannoFlowTypeCombo.currentText())
        
        mValid = bool(self.fannoUpstreamMachEdit.text())
        ttValid = bool(self.fannoTTStEdit.text())
        ppValid = bool(self.fannoPPStEdit.text())
        rrValid = bool(self.fannoRhoRhoStEdit.text())
        popoValid = bool(self.fannoPoPoStEdit.text())
        fricValid = bool(self.fanno4FLStDEdit.text())
        uuValid = bool(self.fannoUUStEdit.text())
        
        state = None
        
        opts = {
            "gamma": gamma,
            "flowType": shockType,
        }
        
        if choice == FannoFlowChoice.GAMMA_MACH and mValid:
            m = float(self.fannoUpstreamMachEdit.text())
            state = FFR(**opts, mach=m)
        
        elif choice == FannoFlowChoice.GAMMA_T_T_ST and ttValid:
            tt = float(self.fannoTTStEdit.text())
            state = FFR(**opts, t_tSt=tt)

        elif choice == FannoFlowChoice.GAMMA_P_P_ST and ppValid:
            pp = float(self.fannoPPStEdit.text())
            state = FFR(**opts, p_pSt=pp)
        
        elif choice == FannoFlowChoice.GAMMA_RHO_RHO_ST and rrValid:
            rr = float(self.fannoRhoRhoStEdit.text())
            state = FFR(**opts, rho_rhoSt=rr)
        
        elif choice == FannoFlowChoice.GAMMA_PO_PO_ST and popoValid:
            popo = float(self.fannoPoPoStEdit.text())
            state = FFR(**opts, po_poSt=popo)
        
        elif choice == FannoFlowChoice.GAMMA_4FLSTD_FLOW_TYPE and fricValid:
            fric = float(self.fanno4FLStDEdit.text())
            state = FFR(**opts, f4LSt_D=fric)
        
        elif choice == FannoFlowChoice.GAMMA_U_U_ST and uuValid:
            uu = float(self.fannoUUStEdit.text())
            state = FFR(**opts, u_uSt=uu)
        
        if state:
            self.fannoUpstreamMachEdit.setText(TO_STR(state.mach))
            self.fannoTTStEdit.setText(TO_STR(state.mach))
            self.fannoPPStEdit.setText(TO_STR(state.p_pSt))
            self.fannoRhoRhoStEdit.setText(TO_STR(state.rho_rhoSt))
            self.fannoPoPoStEdit.setText(TO_STR(state.po_poSt))
            self.fanno4FLStDEdit.setText(TO_STR(state.f4LSt_D))
            self.fannoUUStEdit.setText(TO_STR(state.u_uSt))

        self._fannoState = state

    def calculateFannoFrictionAddition(self) -> None:
        if self._fannoState is None:
            return
        
        validParams = bool(self.fannoPipeDiameterEdit.text())
        validParams &= bool(self.fannoFrictionCoeffEdit.text())
        validParams &= bool(self.fannoPipeLenEdit.text())
        
        if not validParams:
            return
        
        d = float(self.fannoPipeDiameterEdit.text())
        f = float(self.fannoFrictionCoeffEdit.text())
        l = float(self.fannoPipeLenEdit.text())
        self._fannoState.apply_pipe_parameters(d, l, f)
        
        s = self._fannoState
        
        # Down stream conditions
        self.fannoDwnStrmMachEdit.setText(TO_STR(s.dwnStrmMach))
        self.fannoDwnStrmTTStEdit.setText(TO_STR(s.dwnStrm_t_tSt))
        self.fannoDwnStrmPPStEdit.setText(TO_STR(s.dwnStrm_p_pSt))
        self.fannoDwnStrmRhoRhoStEdit.setText(TO_STR(s.dwnStrm_rho_rhoSt))
        self.fannoDwnStrmPoPoStrEdit.setText(TO_STR(s.dwnStrm_f4LSt_D))
        self.fannoDwnStrmUUStEdit.setText(TO_STR(s.dwnStrm_u_uSt))
        self.fannoChokedLengthEdit.setText(TO_STR(s.chokedLength))
        
        # Ratio of down stream to upstream conditions
        self.fannoT2T1Edit.setText(TO_STR(s.t2_t1))
        self.fannoP2P1Edit.setText(TO_STR(s.p2_p1))
        self.fannoU2U1Edit.setText(TO_STR(s.u2_u1))
        self.fannoRho2Rho1Edit.setText(TO_STR(s.rho2_rho1))
        self.fannoPo2Po1Edit.setText(TO_STR(s.po2_po1))
        self.fanno4FLStD24FLStD1Edit.setText(TO_STR(s.f4LD2_f4LD1))
        

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