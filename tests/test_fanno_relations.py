from pytest import approx

from CompAero.fanno_flow_relations import FannoFlowRelations as ffr
from CompAero.internal import FlowState as FS


class TestFannoClassFuncs:
    gamma = 1.4

    # Test the Functions for Subsonic Case
    #######################################################################################
    def test_subsonic_t_tstar(self):
        assert ffr.calc_t_t_star(0.5, self.gamma) == approx(1.1429, rel=1e-4)

    def test_subsonic_mach_from_t_tstar(self):
        assert ffr.calc_mach_from_t_t_star(1.14285714, self.gamma) == approx(0.5, rel=1e-2)

    def test_subsonic_p_pstar(self):
        assert ffr.calc_p_p_star(0.5, self.gamma) == approx(2.1381, rel=1e-4)

    def test_subsonic_mach_from_p_pstar(self):
        assert ffr.calc_mach_from_p_p_star(2.13808993, self.gamma) == approx(0.5, rel=1e-2)

    def test_subsonic_rho_rho_star(self):
        assert ffr.calc_rho_rho_star(0.5, self.gamma) == approx(1.871, rel=1e-4)

    def test_subsonic_mach_from_rho_rho_star(self):
        assert ffr.calc_mach_from_rho_rho_star(1.871, 1.4) == approx(0.5, 1e-3)

    def test_subsonic_p0_p0Star(self):
        assert ffr.calc_po_po_star(0.5, self.gamma) == approx(1.3398, rel=1e-4)

    def test_subsonic_mach_from_p0_p0Star(self):
        assert ffr.calc_mach_from_po_po_star(1.33984375, self.gamma, flow_type=FS.SUB_SONIC) == approx(
            0.5, 1e-3
        )

    def test_subsonic_4FLstarD(self):
        assert ffr.calc_4flst_d(0.5, self.gamma) == approx(1.0691, rel=1e-4)

    def test_subsonic_mach_from_4FLstarD(self):
        assert ffr.calc_mach_from_4flst_d(1.06906031, self.gamma, flow_type=FS.SUB_SONIC) == approx(
            0.5, rel=1e-3
        )

    def test_subsonic_u_u_star(self):
        assert ffr.calc_u_u_starar(0.5, self.gamma) == approx(0.5345, rel=1e-4)

    def test_subsonic_mach_from_u_u_star(self):
        assert ffr.calc_mach_from_u_u_star(0.53452248, self.gamma) == approx(0.5, rel=1e-3)

    # Test the Functions for Supersonic Case
    #######################################################################################
    def test_supersonic_t_tstar(self):
        assert ffr.calc_t_t_star(1.5, self.gamma) == approx(0.82759, rel=1e-4)

    def test_supersonic_mach_from_t_tstar(self):
        assert ffr.calc_mach_from_t_t_star(0.82758620, self.gamma) == approx(1.5, rel=1e-2)

    def test_supersonic_p_pstar(self):
        assert ffr.calc_p_p_star(1.5, self.gamma) == approx(0.6065, rel=1e-4)

    def test_supersonic_mach_from_p_pstar(self):
        assert ffr.calc_mach_from_p_p_star(0.60647843, self.gamma) == approx(1.5, rel=1e-2)

    def test_supersonic_rho_rho_star(self):
        assert ffr.calc_rho_rho_star(1.5, self.gamma) == approx(0.7328, rel=1e-4)

    def test_supersonic_mach_from_rho_rho_star(self):
        assert ffr.calc_mach_from_rho_rho_star(0.7328, 1.4) == approx(1.5, 1e-3)

    def test_supersonic_p0_p0Star(self):
        assert ffr.calc_po_po_star(1.5, self.gamma) == approx(1.1762, rel=1e-4)

    def test_supersonic_mach_from_p0_p0Star(self):
        assert ffr.calc_mach_from_po_po_star(1.17616705, self.gamma, flow_type=FS.SUPER_SONIC) == approx(
            1.5, 1e-3
        )

    def test_supersonic_4FLstarD(self):
        assert ffr.calc_4flst_d(1.5, self.gamma) == approx(0.13605, rel=1e-4)

    def test_supersonic_mach_from_4FLstarD(self):
        assert ffr.calc_mach_from_4flst_d(0.13605021, self.gamma, flow_type=FS.SUPER_SONIC) == approx(
            1.5, rel=1e-3
        )

    def test_supersonic_u_u_star(self):
        assert ffr.calc_u_u_starar(1.5, self.gamma) == approx(1.3646, rel=1e-4)

    def test_supersonic_mach_from_u_u_star(self):
        assert ffr.calc_mach_from_u_u_star(1.36457647, self.gamma) == approx(1.5, rel=1e-3)


class TestFannoClassSubsonic:
    gamma = 1.4

    def test_fanno_from_mach(self):
        inst = ffr(self.gamma, mach=0.5)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_t_t_star(self):
        inst = ffr(self.gamma, t_t_st=1.1428571428571428)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_p_p_star(self):
        inst = ffr(self.gamma, p_p_st=2.1381)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_rho_rho_star(self):
        inst = ffr(self.gamma, rho_rho_st=1.8708286933869707)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_po_po_star(self):
        inst = ffr(self.gamma, po_po_st=1.33984375, flow_type=FS.SUB_SONIC)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_f4LStar_D(self):
        inst = ffr(self.gamma, f4lst_d=1.0690603127182559, flow_type=FS.SUB_SONIC)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_u_u_star(self):
        inst = ffr(self.gamma, u_u_st=0.5345224838248488)
        inst.apply_pipe_parameters(0.4, 11, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(0.593, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.1211, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.7855, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1966, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.5926, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.5191, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.6279, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_choked_flow(self):
        inst = ffr(self.gamma, mach=0.5)
        inst.apply_pipe_parameters(0.4, 22, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(0.5, rel=1e-3)
        assert inst.t_t_st == approx(1.1429, rel=1e-4)
        assert inst.p_p_st == approx(2.1381, rel=1e-4)
        assert inst.rho_rho_st == approx(1.871, rel=1e-4)
        assert inst.po_po_st == approx(1.3398, rel=1e-4)
        assert inst.f4lst_d == approx(1.0691, rel=1e-4)
        assert inst.u_u_st == approx(0.5345, rel=1e-4)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.0, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.0, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.0, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.0, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)


class TestFannoClassSupersonic:
    gamma = 1.4

    def test_fanno_from_mach(self):
        inst = ffr(self.gamma, mach=1.5)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_t_t_star(self):
        inst = ffr(self.gamma, t_t_st=0.8275862068965517)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_p_p_star(self):
        inst = ffr(self.gamma, p_p_st=0.6064784348631227)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_rho_rho_star(self):
        inst = ffr(self.gamma, rho_rho_st=0.7328281087929399)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_po_po_star(self):
        inst = ffr(self.gamma, po_po_st=1.1761670524691357, flow_type=FS.SUPER_SONIC)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_f4LStar_D(self):
        inst = ffr(self.gamma, f4lst_d=0.13605021738414635, flow_type=FS.SUPER_SONIC)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_from_u_u_star(self):
        inst = ffr(self.gamma, u_u_st=1.364576478442026)
        inst.apply_pipe_parameters(0.4, 1.5, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert not inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.2887, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(0.9008, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(0.7365, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0616, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.8176, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.06105, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2231, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)

    def test_fanno_choked_flow(self):
        inst = ffr(self.gamma, mach=1.5)
        inst.apply_pipe_parameters(0.4, 22, 0.005)
        assert inst.gamma == approx(self.gamma, rel=1e-3)
        assert inst.mach == approx(1.5, rel=1e-3)
        assert inst.t_t_st == approx(0.8276, rel=1e-4)
        assert inst.p_p_st == approx(0.6065, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7328, rel=1e-4)
        assert inst.po_po_st == approx(1.1762, rel=1e-4)
        assert inst.f4lst_d == approx(0.13605, rel=1e-4)
        assert inst.u_u_st == approx(1.3646, rel=1e-4)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.choked_flow
        assert inst.dwn_strm_mach == approx(1.0, 1e-3)
        assert inst.dwn_strm_t_t_st == approx(1.0, 1e-3)
        assert inst.dwn_strm_p_p_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_f4lst_d == approx(0.0, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.0, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-5)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-5)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-5)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-5)
        assert inst.f4ld2_f4ld1 == approx(inst.dwn_strm_f4lst_d / inst.f4lst_d, rel=1e-5)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-5)
