"""
This module encompasses the entire compaero ui
"""


import sys
from math import nan
from typing import Optional, TypedDict, Union

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow

from CompAero.Calculator.CalculatorUI import Ui_MainWindow
from CompAero.Calculator.decorators import error_message_decorator
from CompAero.fanno_flow_relations import FANNO_FLOW_VALID_OPTIONS, FannoFlowChoice
from CompAero.fanno_flow_relations import FannoFlowRelations as FFR
from CompAero.internal import FlowState, ShockType
from CompAero.isentropic_relations import (
    ISENTROPIC_CHOICE,
    ISENTROPIC_VALID_OPTIONS,
    IsentropicRelations,
)
from CompAero.normal_shock_relations import NORMAL_SHOCK_CHOICE, NORMAL_SHOCK_VALID_OPTIONS
from CompAero.normal_shock_relations import NormalShockRelations as NSR
from CompAero.oblique_shock_relations import OBLIQUE_SHOCK_VALID_OPTIONS, ObliqueShockChoice
from CompAero.oblique_shock_relations import ObliqueShockRelations as OSR
from CompAero.prandtl_meyer import PRANDTL_MEYER_OPTIONS
from CompAero.prandtl_meyer import PrandtlMeyer as PM
from CompAero.prandtl_meyer import PrandtlMeyerChoice
from CompAero.rayleigh_flow_relations import RAYLEIGH_FLOW_VALID_OPTIONS, RayleighFlowChoice
from CompAero.rayleigh_flow_relations import RayleighFlowRelations as RFR
from CompAero.rocket_nozzle import (
    max_thrust_coefficient,
    min_thrust_coefficient,
    thrust_coefficient,
)

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)  # type: ignore
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)  # type: ignore

PRECISION = 6


def to_str(value: Union[int, float], decimal: int = PRECISION) -> str:
    """
    This funciton rounds the value and converts in to a string

    Args:
        decimal (int): how many digits to round the value to
        value (float | int): value to be rounded and converted to a string

    Returns: a string version of value rounded to decimal digits
    """
    return str(round(value, decimal))


PRANDTL_MEYER_DICT_TYPE = TypedDict(
    "PRANDTL_MEYER_DICT_TYPE", {"gamma": float, "in_degrees": bool, "deflection_angle": float}
)
RAY_FANNO_DICT_TYPE = TypedDict("RAY_FANNO_DICT_TYPE", {"gamma": float, "flow_type": FlowState})
OBLIQUE_SHOCK_REL_DICT_TYPE = TypedDict(
    "OBLIQUE_SHOCK_REL_DICT_TYPE",
    {"gamma": float, "use_degrees": bool, "shock_type": ShockType, "shock_angle": float},
)


class UI(QMainWindow, Ui_MainWindow):
    """Class defines the comp aero ui"""

    def __init__(self) -> None:
        super().__init__()
        self.setupUi(self)
        self.show()
        self.layout().setContentsMargins(0, 0, 0, 0)  # type: ignore

        self._error_box = QtWidgets.QErrorMessage()

        # Combos
        self.isentropicOptionCombo.addItems(ISENTROPIC_VALID_OPTIONS)
        self.isentropicFlowTypeCombo.addItems(
            [FlowState.SUPER_SONIC.name, FlowState.SUB_SONIC.name]
        )
        self.normalShockOptionCombo.addItems(NORMAL_SHOCK_VALID_OPTIONS)
        self.obliqueShockOptionCombo.addItems(OBLIQUE_SHOCK_VALID_OPTIONS)
        self.obliqueshock_typeCombo.addItems([x.value for x in ShockType])
        self.fannoOptionCombo.addItems(FANNO_FLOW_VALID_OPTIONS)
        self.fannoFlowTypeCombo.addItems([x.value for x in FlowState])
        self.rayleighOptionCombo.addItems(RAYLEIGH_FLOW_VALID_OPTIONS)
        self.rayleighFlowTypeCombo.addItems([x.value for x in FlowState])
        self.prandtlMeyerOptionCombo.addItems(PRANDTL_MEYER_OPTIONS)

        # Buttons
        self.isentropicCalcBtn.clicked.connect(self.calc_isentropic_state)
        self.normalShockCalculate.clicked.connect(self.calc_normal_shock_state)
        self.obliqueShockCalcBtn.clicked.connect(self.calc_oblique_shock_state)
        self.obliqueShockDegreeChkBtn.clicked.connect(self.oblique_shock_degrees_chk_box_update)
        self.prandtlMeyerDegreeChkBtn.clicked.connect(self.prandtl_meyer_degrees_chk_box_update)
        self.fannoCalculateBtn.clicked.connect(self.calc_fanno_flow_state)
        self.fannoApplyPipeParamBtn.clicked.connect(self.calc_fanno_friction_addition)
        self.rayleighCalculateBtn.clicked.connect(self.calc_rayleigh_flow_state)
        self.rayleighApplyPipeParamBtn.clicked.connect(self.calc_rayleigh_heat_addition)
        self.prandtlMeyerCalculateBtn.clicked.connect(self.calc_prantl_meyer_state)
        self.rocketNozzleCalculateBtn.clicked.connect(self.calc_rocket_nozzle_state)

        # Saved states
        self._fanno_state: Optional[FFR] = None
        self._rayleigh_state: Optional[RFR] = None

    def show_error(self, msg: str) -> None:
        """
        Shows the error popup with a message

        Args:
            msg (str): msg to show in the popup
        """
        self._error_box.showMessage(msg)

    def oblique_shock_degrees_chk_box_update(self) -> None:
        """Changes the text on the button to reflect the oblique shock angle units"""
        if self.obliqueShockDegreeChkBtn.isChecked():
            self.obliqueShockDegreeChkBtn.setText("Radians")
        else:
            self.obliqueShockDegreeChkBtn.setText("Degrees")

    @error_message_decorator
    def prandtl_meyer_degrees_chk_box_update(self) -> None:
        """Changes the text on the button to reflect the prantl meyer angle units"""
        if self.prandtlMeyerDegreeChkBtn.isChecked():
            self.prandtlMeyerDegreeChkBtn.setText("Radians")
        else:
            self.prandtlMeyerDegreeChkBtn.setText("Degrees")

    @error_message_decorator
    def calc_isentropic_state(self) -> None:
        """Calculates the isentropic state and sets the entry box values"""
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
            choice == ISENTROPIC_CHOICE.a_a_star
            and self.isentropicAAStarEntry.text()
            and self.isentropicFlowTypeCombo.currentText()
        ):
            flowtype = FlowState(self.isentropicFlowTypeCombo.currentText())
            aa_star = float(self.isentropicAAStarEntry.text())
            state = IsentropicRelations(gamma, a_a_star=aa_star, flow_type=flowtype)

        if state is not None:
            self.isentropicMachEntry.setText(to_str(state.mach))
            self.isentropicP0PEntry.setText(to_str(state.p0_p))
            self.isentropicAAStarEntry.setText(to_str(state.a_a_star))
            self.isentropicT0TEntry.setText(to_str(state.t0_t))
            self.isentropicRho0RhoEntry.setText(to_str(state.rho0_rho))

    @error_message_decorator
    def calc_normal_shock_state(self) -> None:
        """Calculates the normal shock state and sets the entry box values"""
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
            self.normalShockM1Entry.setText(to_str(state.mach))
            self.normalShockM2Entry.setText(to_str(state.mach2))
            self.normalShockP2P1Entry.setText(to_str(state.p2_p1))
            self.normalShockRho2Rho1Entry.setText(to_str(state.rho2_rho1))
            self.normalShockT2T1Entry.setText(to_str(state.t2_t1))
            self.normalShockP02P01Entry.setText(to_str(state.po2_po1))
            self.normalShockP02P1Entry.setText(to_str(state.po2_p1))

    @error_message_decorator
    def calc_oblique_shock_state(self) -> None:
        # pylint: disable=too-many-statements
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-locals
        """Calculates the oblique shock state and sets the entry box values"""
        if not self.obliqueShockGammaEntry.text():
            return

        gamma = float(self.obliqueShockGammaEntry.text())
        choice = ObliqueShockChoice(self.obliqueShockOptionCombo.currentText())
        shock_type = ShockType(self.obliqueshock_typeCombo.currentText())
        wedge_angle_valid = bool(self.obliqueShockwedge_angleEdit.text())
        shock_angle_valid = bool(self.obliqueshock_angleEdit.text())
        m1_valid = bool(self.obliqueShockM1Edit.text())
        mn1_valid = bool(self.obliqueShockMn1Edit.text())
        mn2_valid = bool(self.obiqueShockMn2Edit.text())
        m2_valid = bool(self.obliqueShockM2Edit.text())
        p2p1_valid = bool(self.obliqueShockP2P1Edit.text())
        rho2_rho1_valid = bool(self.obliqueShockRho2Rho1Edit.text())
        po2_p1_valid = bool(self.obliqueShockPo2P1Edit.text())
        t2_t1_valid = bool(self.obliqueShockT2T1Edit.text())
        po2_po1_valid = bool(self.obliqueShockPo2Po1Edit.text())
        degrees = bool(not self.obliqueShockDegreeChkBtn.isChecked())

        if not choice:
            return

        state = None

        # common options
        opts: OBLIQUE_SHOCK_REL_DICT_TYPE = {
            "gamma": gamma,
            "use_degrees": degrees,
            "shock_type": shock_type,
            "shock_angle": nan,
        }

        if choice == ObliqueShockChoice.MACH_WEDGE_ANGLE and m1_valid and wedge_angle_valid:
            m = float(self.obliqueShockM1Edit.text())
            wa = float(self.obliqueShockwedge_angleEdit.text())
            state = OSR(**opts, mach=m, wedge_angle=wa)

        if not shock_angle_valid:
            return

        opts["shock_angle"] = float(self.obliqueshock_angleEdit.text())

        if choice == ObliqueShockChoice.MACH_SHOCK_ANGLE and m1_valid:
            m = float(self.obliqueShockM1Edit.text())
            state = OSR(**opts, mach=m)

        elif choice == ObliqueShockChoice.MACH_N_1_SHOCK_ANGLE and mn1_valid:
            m = float(self.obliqueShockMn1Edit.text())
            state = OSR(**opts, mn1=m)

        elif choice == ObliqueShockChoice.M2_WEDGE_SHOCK_ANGLE and m2_valid and wedge_angle_valid:
            m2 = float(self.obliqueShockM2Edit.text())
            wa = float(self.obliqueShockwedge_angleEdit.text())
            state = OSR(**opts, m2=m2, wedge_angle=wa)

        elif choice == ObliqueShockChoice.P2_P1_SHOCK_ANGLE and p2p1_valid:
            p2p1 = float(self.obliqueShockP2P1Edit.text())
            state = OSR(**opts, p2_p1=p2p1)

        elif choice == ObliqueShockChoice.RHO2_RHO1_SHOCK_ANGLE and rho2_rho1_valid:
            r2r1 = float(self.obliqueShockRho2Rho1Edit.text())
            state = OSR(**opts, rho2_rho1=r2r1)

        elif choice == ObliqueShockChoice.T2_T1_SHOCK_ANGLE and t2_t1_valid:
            t2t1 = float(self.obliqueShockT2T1Edit.text())
            state = OSR(**opts, t2_t1=t2t1)

        elif choice == ObliqueShockChoice.PO2_PO1_SHOCK_ANGLE and po2_po1_valid:
            p = float(self.obliqueShockPo2Po1Edit.text())
            state = OSR(**opts, po2_po1=p)

        elif choice == ObliqueShockChoice.PO2_P1_SHOCK_ANGLE and po2_p1_valid:
            p = float(self.obliqueShockPo2P1Edit.text())
            state = OSR(**opts, po2_p1=p)

        elif choice == ObliqueShockChoice.MN2_SHOCK_ANGLE and mn2_valid:
            m = float(self.obiqueShockMn2Edit.text())
            state = OSR(**opts, mn2=m)

        if state:
            self.obliqueShockwedge_angleEdit.setText(to_str(state.wedge_angle))
            self.obliqueshock_angleEdit.setText(to_str(state.shock_angle))
            self.obliqueShockM1Edit.setText(to_str(state.mach))
            self.obliqueShockMn1Edit.setText(to_str(state.mach_normal_1))
            self.obiqueShockMn2Edit.setText(to_str(state.mach_normal_2))
            self.obliqueShockM2Edit.setText(to_str(state.mach2))
            self.obliqueShockP2P1Edit.setText(to_str(state.p2_p1))
            self.obliqueShockRho2Rho1Edit.setText(to_str(state.rho2_rho1))
            self.obliqueShockPo2P1Edit.setText(to_str(state.po2_p1))
            self.obliqueShockT2T1Edit.setText(to_str(state.t2_t1))
            self.obliqueShockPo2Po1Edit.setText(to_str(state.po2_po1))

    @error_message_decorator
    def calc_fanno_flow_state(self) -> None:
        # pylint: disable=too-many-locals
        """Calculates the fanno flow state and sets the entry box values"""
        if not self.fannoGammaEntry.text():
            return

        gamma = float(self.fannoGammaEntry.text())
        choice = FannoFlowChoice(self.fannoOptionCombo.currentText())
        flow_type: FlowState = FlowState(self.fannoFlowTypeCombo.currentText())

        m_valid = bool(self.fannoUpstreamMachEdit.text())
        tt_valid = bool(self.fannoTTStEdit.text())
        pp_valid = bool(self.fannoPPStEdit.text())
        rr_valid = bool(self.fannoRhoRhoStEdit.text())
        popo_valid = bool(self.fannoPoPoStEdit.text())
        fric_valid = bool(self.fanno4FLStDEdit.text())
        uu_valid = bool(self.fannoUUStEdit.text())

        state = None

        opts: RAY_FANNO_DICT_TYPE = {"gamma": gamma, "flow_type": flow_type}

        if choice == FannoFlowChoice.GAMMA_MACH and m_valid:
            m = float(self.fannoUpstreamMachEdit.text())
            state = FFR(**opts, mach=m)

        elif choice == FannoFlowChoice.GAMMA_T_T_ST and tt_valid:
            tt = float(self.fannoTTStEdit.text())
            state = FFR(**opts, t_t_st=tt)

        elif choice == FannoFlowChoice.GAMMA_P_P_ST and pp_valid:
            pp = float(self.fannoPPStEdit.text())
            state = FFR(**opts, p_p_st=pp)

        elif choice == FannoFlowChoice.GAMMA_RHO_RHO_ST and rr_valid:
            rr = float(self.fannoRhoRhoStEdit.text())
            state = FFR(**opts, rho_rho_st=rr)

        elif choice == FannoFlowChoice.GAMMA_PO_PO_ST and popo_valid:
            popo = float(self.fannoPoPoStEdit.text())
            state = FFR(**opts, po_po_st=popo)

        elif choice == FannoFlowChoice.GAMMA_4FLSTD_FLOW_TYPE and fric_valid:
            fric = float(self.fanno4FLStDEdit.text())
            state = FFR(**opts, f4lst_d=fric)

        elif choice == FannoFlowChoice.GAMMA_U_U_ST and uu_valid:
            uu = float(self.fannoUUStEdit.text())
            state = FFR(**opts, u_u_st=uu)

        if state:
            self.fannoUpstreamMachEdit.setText(to_str(state.mach))
            self.fannoTTStEdit.setText(to_str(state.mach))
            self.fannoPPStEdit.setText(to_str(state.p_p_st))
            self.fannoRhoRhoStEdit.setText(to_str(state.rho_rho_st))
            self.fannoPoPoStEdit.setText(to_str(state.po_po_st))
            self.fanno4FLStDEdit.setText(to_str(state.f4lst_d))
            self.fannoUUStEdit.setText(to_str(state.u_u_st))

        self._fanno_state = state

    @error_message_decorator
    def calc_fanno_friction_addition(self) -> None:
        """Calculates the addition of friction to the fanno state and sets the entry box values"""
        if self._fanno_state is None:
            return

        valid_params = bool(self.fannoPipeDiameterEdit.text())
        valid_params &= bool(self.fannoFrictionCoeffEdit.text())
        valid_params &= bool(self.fannoPipeLenEdit.text())

        if not valid_params:
            return

        d = float(self.fannoPipeDiameterEdit.text())
        f = float(self.fannoFrictionCoeffEdit.text())
        l = float(self.fannoPipeLenEdit.text())
        self._fanno_state.apply_pipe_parameters(d, l, f)

        s = self._fanno_state

        # Down stream conditions
        self.fannoDwnStrmMachEdit.setText(to_str(s.dwn_strm_mach))
        self.fannoDwnStrmTTStEdit.setText(to_str(s.dwn_strm_t_t_st))
        self.fannoDwnStrmPPStEdit.setText(to_str(s.dwn_strm_p_p_st))
        self.fannoDwnStrmRhoRhoStEdit.setText(to_str(s.dwn_strm_rho_rho_st))
        self.fannoDwnStrmPoPoStrEdit.setText(to_str(s.dwn_strm_f4lst_d))
        self.fannoDwnStrmUUStEdit.setText(to_str(s.dwn_strm_u_u_st))
        self.fannoChokedLengthEdit.setText(to_str(s.choked_length))

        # Ratio of down stream to upstream conditions
        self.fannoT2T1Edit.setText(to_str(s.t2_t1))
        self.fannoP2P1Edit.setText(to_str(s.p2_p1))
        self.fannoU2U1Edit.setText(to_str(s.u2_u1))
        self.fannoRho2Rho1Edit.setText(to_str(s.rho2_rho1))
        self.fannoPo2Po1Edit.setText(to_str(s.po2_po1))
        self.fanno4FLStD24FLStD1Edit.setText(to_str(s.f4ld2_f4ld1))

    @error_message_decorator
    def calc_rayleigh_flow_state(self) -> None:
        """Calculates the rayleigh flow state and sets the entry box values"""
        if not self.rayleighGammaEntry.text():
            return

        gamma = float(self.rayleighGammaEntry.text())
        flow_type = FlowState(self.rayleighFlowTypeCombo.currentText())
        choice = RayleighFlowChoice(self.rayleighOptionCombo.currentText())

        m_valid = bool(self.rayleighUpstreamMachEdit.text())
        tt_valid = bool(self.rayleighTTStEdit.text())
        pp_valid = bool(self.rayleighPPStEdit.text())
        rr_valid = bool(self.rayleighRhoRhoStEdit.text())
        popo_valid = bool(self.rayleighPoPoStEdit.text())
        toto_valid = bool(self.rayleighToToStEdit.text())
        uu_valid = bool(self.rayleighUUStEdit.text())

        state = None

        opts: RAY_FANNO_DICT_TYPE = {"gamma": gamma, "flow_type": flow_type}

        if choice == RayleighFlowChoice.GAMMA_MACH and m_valid:
            state = RFR(**opts, mach=float(self.rayleighUpstreamMachEdit.text()))

        elif choice == RayleighFlowChoice.GAMMA_T_T_ST and tt_valid:
            state = RFR(**opts, t_t_st=float(self.rayleighTTStEdit.text()))

        elif choice == RayleighFlowChoice.GAMMA_P_P_ST and pp_valid:
            state = RFR(**opts, p_p_st=float(self.rayleighPPStEdit.text()))

        elif choice == RayleighFlowChoice.GAMMA_RHO_RHO_ST and rr_valid:
            state = RFR(**opts, rho_rho_st=float(self.rayleighRhoRhoStEdit.text()))

        elif choice == RayleighFlowChoice.GAMMA_PO_PO_ST and popo_valid:
            state = RFR(**opts, po_po_st=float(self.rayleighPoPoStEdit.text()))

        elif choice == RayleighFlowChoice.GAMMA_TO_TO_FLOW_TYPE and toto_valid:
            state = RFR(**opts, to_to_st=float(self.rayleighToToStEdit.text()))

        elif choice == RayleighFlowChoice.GAMMA_U_U_ST and uu_valid:
            state = RFR(**opts, u_u_st=float(self.rayleighUUStEdit.text()))

        self._rayleigh_state = state

        if state:
            s = state
            self.rayleighUpstreamMachEdit.setText(to_str(s.mach))
            self.rayleighTTStEdit.setText(to_str(s.t_t_st))
            self.rayleighPPStEdit.setText(to_str(s.p_p_st))
            self.rayleighRhoRhoStEdit.setText(to_str(s.rho_rho_st))
            self.rayleighPoPoStEdit.setText(to_str(s.po_po_st))
            self.rayleighToToStEdit.setText(to_str(s.to_to_st))
            self.rayleighUUStEdit.setText(to_str(s.u_u_st))

    @error_message_decorator
    def calc_rayleigh_heat_addition(self) -> None:
        """Calculates the addition of heat to the rayleigh state and sets the entry box values"""
        if self._rayleigh_state is None:
            return

        valid_params = bool(self.rayleighHeatEdit.text())
        valid_params &= bool(self.rayleighGasConstantEdit.text())
        valid_params &= bool(self.rayleighHeatTo1Edit.text())

        if not valid_params:
            return

        h = float(self.rayleighHeatEdit.text())
        r = float(self.rayleighGasConstantEdit.text())
        t = float(self.rayleighHeatTo1Edit.text())

        self._rayleigh_state.simulate_heat_addition(h, t, r)
        s = self._rayleigh_state

        # Down stream conditions
        self.rayleighDwnStrmMachEdit.setText(to_str(s.dwn_strm_mach))
        self.rayleighDwnStrmTTStEdit.setText(to_str(s.dwn_strm_t_t_st))
        self.rayleighDwnStrmPPStEdit.setText(to_str(s.dwn_strm_p_p_st))
        self.rayleighDwnStrmRhoRhoStEdit.setText(to_str(s.dwn_strm_rho_rho_st))
        self.rayleighDwnStrmPoPoStrEdit.setText(to_str(s.dwn_strm_po_po_st))
        self.rayleighDwnStrmToToStEdit.setText(to_str(s.dwn_strm_to_to_st))
        self.rayleighDwnStrmUUStEdit.setText(to_str(s.dwn_strm_u_u_st))
        self.rayleighChokedHeatEdit.setText(to_str(s.choked_heat))

        # Ratio of down stream to upstream conditions
        self.rayleighT2T1Edit.setText(to_str(s.t2_t1))
        self.rayleighP2P1Edit.setText(to_str(s.p2_p1))
        self.rayleighU2U1Edit.setText(to_str(s.u2_u1))
        self.rayleighRho2Rho1Edit.setText(to_str(s.rho2_rho1))
        self.rayleighPo2Po1Edit.setText(to_str(s.po2_po1))
        self.rayleighTo2To1Edit.setText(to_str(s.to2_to1))
        self.rayleighTo2Edit.setText(to_str(s.to2))

    @error_message_decorator
    def calc_prantl_meyer_state(self) -> None:
        """Calculates prantl meyer state and sets the entry box values"""
        if not self.prandtlMeyerGammaEntry.text():
            return

        gamma = float(self.prandtlMeyerGammaEntry.text())

        choice = PrandtlMeyerChoice(self.prandtlMeyerOptionCombo.currentText())

        m_valid = bool(self.prandtlMeyerUpstreamMachEdit.text())
        nu_valid = bool(self.prandtlMeyerNuEdit.text())
        mu_valid = bool(self.prandtlMeyerMuEdit.text())
        angle_valid = bool(self.prandtlMeyerDeflectionAngleEdit.text())
        d_m_valid = bool(self.prandtlMeyerDwnStrmMachEdit.text())
        d_nu_valid = bool(self.prandtlMeyerDwnStrmNuEdit.text())
        d_mu_valid = bool(self.prandtlMeyerDwnStrmMuEdit.text())
        degrees = not self.prandtlMeyerDegreeChkBtn.isChecked()

        state = None

        opts: PRANDTL_MEYER_DICT_TYPE = {
            "gamma": gamma,
            "in_degrees": degrees,
            "deflection_angle": nan,
        }

        if angle_valid:
            opts["deflection_angle"] = float(self.prandtlMeyerDeflectionAngleEdit.text())

        if choice == PrandtlMeyerChoice.GAMMA_MACH and m_valid:
            state = PM(**opts, mach=float(self.prandtlMeyerUpstreamMachEdit.text()))

        elif choice == PrandtlMeyerChoice.GAMMA_NU and nu_valid:
            state = PM(**opts, nu=float(self.prandtlMeyerNuEdit.text()))

        elif choice == PrandtlMeyerChoice.GAMMA_MU and mu_valid:
            state = PM(**opts, mu=float(self.prandtlMeyerMuEdit.text()))

        elif choice == PrandtlMeyerChoice.GAMMA_DEFLECTION_DWN_STRM_MACH and d_m_valid:
            state = PM(**opts, down_stream_mach=float(self.prandtlMeyerDwnStrmMachEdit.text()))

        elif choice == PrandtlMeyerChoice.GAMMA_DEFLECTION_DWN_STRM_MU and d_mu_valid:
            state = PM(**opts, down_stream_mu=float(self.prandtlMeyerDwnStrmMuEdit.text()))

        elif choice == PrandtlMeyerChoice.GAMMA_DEFLECTION_DWN_STRM_NU and d_nu_valid:
            state = PM(**opts, down_stream_nu=float(self.prandtlMeyerDwnStrmNuEdit.text()))

        if state:
            s = state
            self.prandtlMeyerUpstreamMachEdit.setText(to_str(s.mach))
            self.prandtlMeyerNuEdit.setText(to_str(s.nu))
            self.prandtlMeyerMuEdit.setText(to_str(s.mu))
            self.prandtlMeyerDwnStrmMachEdit.setText(to_str(s.down_stream_mach))
            self.prandtlMeyerDwnStrmMuEdit.setText(to_str(s.down_stream_mu))
            self.prandtlMeyerDwnStrmNuEdit.setText(to_str(s.down_stream_nu))

    @error_message_decorator
    def calc_rocket_nozzle_state(self) -> None:
        """Calculates rocket nozzle state and sets the entry box values"""
        gamma_valid = bool(self.rocketNozzleGammaEntry.text())
        ar_valid = bool(self.rocketNozzleAreaRatioEdit.text())
        pepc_valid = bool(self.rocketNozzlePePcEdit.text())
        papc_valid = bool(self.rocketNozzlePaPcEdit.text())

        gamma, ar, pe_pc, pa_pc = (0.0, 0.0, 0.0, 0.0)

        if gamma_valid:
            gamma = float(self.rocketNozzleGammaEntry.text())

        if ar_valid:
            ar = float(self.rocketNozzleAreaRatioEdit.text())

        if pepc_valid:
            pe_pc = float(self.rocketNozzlePePcEdit.text())

        if papc_valid:
            pa_pc = float(self.rocketNozzlePaPcEdit.text())

        if ar_valid:
            min_coeff = min_thrust_coefficient(ar)
            self.rocketNozzleMinThrustCoeffEdit.setText(to_str(min_coeff))

        if gamma_valid and ar_valid and pepc_valid:
            max_coeff = max_thrust_coefficient(gamma, ar, pe_pc)
            self.rocketNozzleMaxThrustCoeffEdit.setText(to_str(max_coeff))

        if gamma_valid and ar_valid and pepc_valid and papc_valid:
            coeff = thrust_coefficient(gamma, ar, pe_pc, pa_pc)
            self.rocketNozzleThrustCoeffEdit.setText(to_str(coeff))


def main() -> None:
    """Main function of application"""
    app = QtWidgets.QApplication(sys.argv)
    _ = UI()
    app.exec_()


if __name__ == "__main__":
    main()
