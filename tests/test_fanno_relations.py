from pytest import approx
from CompAero.FannoFlowRelations import FannoFlowRelations as ffr
from CompAero.internal import FlowState as FS


class TestFannoClassFuncs:
    gamma = 1.4

    # Test the Functions for Subsonic Case
    #######################################################################################
    def test_subsonic_t_tstar(self):
        assert ffr.calc_T_Tstar(0.5, self.gamma) == approx(1.1429, rel=1e-4)

    def test_subsonic_mach_from_t_tstar(self):
        assert ffr.calc_mach_from_T_TStar(1.14285714, self.gamma) == approx(0.5, rel=1e-2)

    def test_subsonic_p_pstar(self):
        assert ffr.calc_P_Pstar(0.5, self.gamma) == approx(2.1381, rel=1e-4)

    def test_subsonic_mach_from_p_pstar(self):
        assert ffr.calc_mach_from_P_PStar(2.13808993, self.gamma) == approx(0.5, rel=1e-2)

    def test_subsonic_rho_rhoStar(self):
        assert ffr.calc_Rho_RhoStar(0.5, self.gamma) == approx(1.871, rel=1e-4)

    def test_subsonic_mach_from_rho_rhoStar(self):
        assert ffr.calc_mach_from_Rho_RhoStar(1.871, 1.4) == approx(0.5, 1e-3)

    def test_subsonic_p0_p0Star(self):
        assert ffr.calc_Po_PoStar(0.5, self.gamma) == approx(1.3398, rel=1e-4)

    def test_subsonic_mach_from_p0_p0Star(self):
        assert ffr.calc_mach_from_Po_PoStar(1.33984375, self.gamma, flowType=FS.SUB_SONIC) == approx(
            0.5, 1e-3
        )

    def test_subsonic_4FLstarD(self):
        assert ffr.calc_4FLSt_D(0.5, self.gamma) == approx(1.0691, rel=1e-4)

    def test_subsonic_mach_from_4FLstarD(self):
        assert ffr.calc_mach_from_4FLSt_D(1.06906031, self.gamma, flowType=FS.SUB_SONIC) == approx(
            0.5, rel=1e-3
        )

    def test_subsonic_u_uStar(self):
        assert ffr.calc_U_UStar(0.5, self.gamma) == approx(0.5345, rel=1e-4)

    def test_subsonic_mach_from_u_uStar(self):
        assert ffr.calc_mach_from_U_USt(0.53452248, self.gamma) == approx(0.5, rel=1e-3)

    # Test the Functions for Supersonic Case
    #######################################################################################
    def test_supersonic_t_tstar(self):
        assert ffr.calc_T_Tstar(1.5, self.gamma) == approx(0.82759, rel=1e-4)

    def test_supersonic_mach_from_t_tstar(self):
        assert ffr.calc_mach_from_T_TStar(0.82758620, self.gamma) == approx(1.5, rel=1e-2)

    def test_supersonic_p_pstar(self):
        assert ffr.calc_P_Pstar(1.5, self.gamma) == approx(0.6065, rel=1e-4)

    def test_supersonic_mach_from_p_pstar(self):
        assert ffr.calc_mach_from_P_PStar(0.60647843, self.gamma) == approx(1.5, rel=1e-2)

    def test_supersonic_rho_rhoStar(self):
        assert ffr.calc_Rho_RhoStar(1.5, self.gamma) == approx(0.7328, rel=1e-4)

    def test_supersonic_mach_from_rho_rhoStar(self):
        assert ffr.calc_mach_from_Rho_RhoStar(0.7328, 1.4) == approx(1.5, 1e-3)

    def test_supersonic_p0_p0Star(self):
        assert ffr.calc_Po_PoStar(1.5, self.gamma) == approx(1.1762, rel=1e-4)

    def test_supersonic_mach_from_p0_p0Star(self):
        assert ffr.calc_mach_from_Po_PoStar(1.17616705, self.gamma, flowType=FS.SUPER_SONIC) == approx(
            1.5, 1e-3
        )

    def test_supersonic_4FLstarD(self):
        assert ffr.calc_4FLSt_D(1.5, self.gamma) == approx(0.13605, rel=1e-4)

    def test_supersonic_mach_from_4FLstarD(self):
        assert ffr.calc_mach_from_4FLSt_D(0.13605021, self.gamma, flowType=FS.SUPER_SONIC) == approx(
            1.5, rel=1e-3
        )

    def test_supersonic_u_uStar(self):
        assert ffr.calc_U_UStar(1.5, self.gamma) == approx(1.3646, rel=1e-4)

    def test_supersonic_mach_from_u_uStar(self):
        assert ffr.calc_mach_from_U_USt(1.36457647, self.gamma) == approx(1.5, rel=1e-3)


class TestFannoClassSubsonic:
    gamma = 1.4

    def test_fanno_from_mach(self):
        inst = ffr(self.gamma, mach=0.5)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_t_tStar(self):
        inst = ffr(self.gamma, t_tSt=1.1428571428571428)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_p_pStar(self):
        inst = ffr(self.gamma, p_pSt=2.1381)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_rho_rhoStar(self):
        inst = ffr(self.gamma, rho_rhoSt=1.8708286933869707)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_po_poStar(self):
        inst = ffr(self.gamma, po_poSt=1.33984375, flowType=FS.SUB_SONIC)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_f4LStar_D(self):
        inst = ffr(self.gamma, f4LSt_D=1.0690603127182559, flowType=FS.SUB_SONIC)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_u_uStar(self):
        inst = ffr(self.gamma, u_uSt=0.5345224838248488)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(0.593, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.1211, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.7855, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1966, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.5926, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.5191, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_choked_flow(self):
        inst = ffr(self.gamma, mach=0.5)
        inst.apply_pipe_parameters(0.4, 22, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_tSt == approx(1.1429, rel=1e-4)
        assert inst.p_pSt == approx(2.1381, rel=1e-4)
        assert inst.rho_rhoSt == approx(1.871, rel=1e-4)
        assert inst.po_poSt == approx(1.3398, rel=1e-4)
        assert inst.f4LSt_D == approx(1.0691, rel=1e-4)
        assert inst.u_uSt == approx(0.5345, rel=1e-4)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.0, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.0, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.0, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.0, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)


class TestFannoClassSupersonic:
    gamma = 1.4

    def test_fanno_from_mach(self):
        inst = ffr(self.gamma, mach=1.5)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_t_tStar(self):
        inst = ffr(self.gamma, t_tSt=0.8275862068965517)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_p_pStar(self):
        inst = ffr(self.gamma, p_pSt=0.6064784348631227)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_rho_rhoStar(self):
        inst = ffr(self.gamma, rho_rhoSt=0.7328281087929399)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_po_poStar(self):
        inst = ffr(self.gamma, po_poSt=1.1761670524691357, flowType=FS.SUPER_SONIC)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_f4LStar_D(self):
        inst = ffr(self.gamma, f4LSt_D=0.13605021738414635, flowType=FS.SUPER_SONIC)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_from_u_uStar(self):
        inst = ffr(self.gamma, u_uSt=1.364576478442026)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert not inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.2887, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(0.9008, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(0.7365, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0616, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.8176, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.06105, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)

    def test_fanno_choked_flow(self):
        inst = ffr(self.gamma, mach=1.5)
        inst.apply_pipe_parameters(0.4, 22, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_tSt == approx(0.8276, rel=1e-4)
        assert inst.p_pSt == approx(0.6065, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7328, rel=1e-4)
        assert inst.po_poSt == approx(1.1762, rel=1e-4)
        assert inst.f4LSt_D == approx(0.13605, rel=1e-4)
        assert inst.u_uSt == approx(1.3646, rel=1e-4)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.chockedFlow
        assert inst.dwnStrmMach == approx(1.0, 1e-3)
        assert inst.dwnStrm_t_tSt == approx(1.0, 1e-3)
        assert inst.dwnStrm_p_pSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_f4LSt_D == approx(0.0, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.0, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-5)
        assert inst.f4LD2_f4LD1 == approx(inst.dwnStrm_f4LSt_D / inst.f4LSt_D, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-5)
