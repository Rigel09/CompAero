# pylint: disable=missing-module-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from pytest import approx

from CompAero.internal import FlowState as FS
from CompAero.rayleigh_flow_relations import RayleighFlowRelations as rfr


class TestRayleighClassFuncs:
    # pylint: disable=too-many-public-methods
    gamma = 1.4

    # Test the Functions for Subsonic Case
    #######################################################################################
    def test_subsonic_t_tstar(self):
        assert rfr.calc_t_t_star(0.5, self.gamma) == approx(0.7901, rel=1e-4)

    def test_subsonic_mach_from_t_tstar(self):
        assert rfr.calc_mach_from_t_t_star(0.79012345, self.gamma, FS.SUB_SONIC) == approx(
            0.5, rel=1e-2
        )

    def test_subsonic_p_pstar(self):
        assert rfr.calc_p_p_star(0.5, self.gamma) == approx(1.7778, rel=1e-4)

    def test_subsonic_mach_from_p_pstar(self):
        assert rfr.calc_mach_from_p_p_star(1.77777777, self.gamma) == approx(0.5, rel=1e-2)

    def test_subsonic_rho_rho_star(self):
        assert rfr.calc_rho_rho_star(0.5, self.gamma) == approx(2.25, rel=1e-4)

    def test_subsonic_mach_from_rho_rho_star(self):
        assert rfr.calc_mach_from_rho_rho_star(2.25, 1.4) == approx(0.5, 1e-3)

    def test_subsonic_p0_p0_star(self):
        assert rfr.calc_po_po_star(0.5, self.gamma) == approx(1.1141, rel=1e-4)

    def test_subsonic_mach_from_p0_p0_star(self):
        assert rfr.calc_mach_from_po_po_star(
            1.11405250, self.gamma, flow_type=FS.SUB_SONIC
        ) == approx(0.5, 1e-3)

    def test_subsonic_t0_t0_star(self):
        assert rfr.calc_to_to_star(0.5, self.gamma) == approx(0.6914, rel=1e-4)

    def test_subsonic_mach_from_t0_t0_star(self):
        assert rfr.calc_mach_from_to_to_star(
            0.69135802, self.gamma, flow_type=FS.SUB_SONIC
        ) == approx(0.5, rel=1e-3)

    def test_subsonic_u_u_star(self):
        assert rfr.calc_u_u_starar(0.5, self.gamma) == approx(0.44444, rel=1e-4)

    def test_subsonic_mach_from_u_u_star(self):
        assert rfr.calc_mach_from_u_u_star(0.44444444, self.gamma) == approx(0.5, rel=1e-3)

    # Test the Functions for Supersonic Case
    #######################################################################################
    def test_supersonic_t_tstar(self):
        assert rfr.calc_t_t_star(1.5, self.gamma) == approx(0.7525, rel=1e-4)

    def test_supersonic_mach_from_t_tstar(self):
        assert rfr.calc_mach_from_t_t_star(
            0.75250399, self.gamma, flow_type=FS.SUPER_SONIC
        ) == approx(1.5, rel=1e-2)

    def test_supersonic_p_pstar(self):
        assert rfr.calc_p_p_star(1.5, self.gamma) == approx(0.5783, rel=1e-4)

    def test_supersonic_mach_from_p_pstar(self):
        assert rfr.calc_mach_from_p_p_star(0.57831325, self.gamma) == approx(1.5, rel=1e-2)

    def test_supersonic_rho_rho_star(self):
        assert rfr.calc_rho_rho_star(1.5, self.gamma) == approx(0.7685, rel=1e-4)

    def test_supersonic_mach_from_rho_rho_star(self):
        assert rfr.calc_mach_from_rho_rho_star(0.7685185185185186, 1.4) == approx(1.5, 1e-3)

    def test_supersonic_p0_p0_star(self):
        assert rfr.calc_po_po_star(1.5, self.gamma) == approx(1.1215, rel=1e-4)

    def test_supersonic_mach_from_p0_p0_star(self):
        assert rfr.calc_mach_from_po_po_star(1.12154522, self.gamma) == approx(1.5, 1e-3)

    def test_supersonic_t0_t0_star(self):
        assert rfr.calc_to_to_star(1.5, self.gamma) == approx(0.9093, rel=1e-4)

    def test_supersonic_mach_from_t0_t0_star(self):
        assert rfr.calc_mach_from_to_to_star(
            0.90927565, self.gamma, flow_type=FS.SUPER_SONIC
        ) == approx(1.5, rel=1e-3)

    def test_supersonic_u_u_star(self):
        assert rfr.calc_u_u_starar(1.5, self.gamma) == approx(1.3012, rel=1e-4)

    def test_supersonic_mach_from_u_u_star(self):
        assert rfr.calc_mach_from_u_u_star(1.30120481, self.gamma) == approx(1.5, rel=1e-3)


class TestRayleighClassSubsonic:
    gamma = 1.4

    def test_rayleigh_from_mach(self):
        inst = rfr(self.gamma, 0.5)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_t_t_star(self):
        inst = rfr(self.gamma, t_t_st=0.7901234567901234, flow_type=FS.SUB_SONIC)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_p_p_star(self):
        inst = rfr(self.gamma, p_p_st=1.7777777777777777)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_rho_rho_star(self):
        inst = rfr(self.gamma, rho_rho_st=2.250)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_po_po_star(self):
        inst = rfr(self.gamma, po_po_st=1.114052503180089, flow_type=FS.SUB_SONIC)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_to_to_star(self):
        inst = rfr(self.gamma, to_to_st=0.691358024691358, flow_type=FS.SUB_SONIC)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_u_u_star(self):
        inst = rfr(self.gamma, u_u_st=0.4444444444444444)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(0.7525, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0148, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.3388, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.3192, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0295, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(0.7580, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_choked_flow(self):
        inst = rfr(self.gamma, u_u_st=0.4444444444444444)
        inst.simulate_heat_addition(200000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flow_type == FS.SUB_SONIC
        assert inst.t_t_st == approx(0.7901, rel=1e-4)
        assert inst.p_p_st == approx(1.7778, rel=1e-4)
        assert inst.rho_rho_st == approx(2.25, rel=1e-2)
        assert inst.po_po_st == approx(1.1141, rel=1e-4)
        assert inst.to_to_st == approx(0.6914, rel=1e-4)
        assert inst.u_u_st == approx(0.44444, rel=1e-4)
        assert inst.choked_heat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(200000)
        assert inst.gas_constant_r == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(1.0, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(474.304, rel=1e-3)
        assert inst.choked_flow


class TestRayleighClassSuperSonic:
    gamma = 1.4

    def test_rayleigh_from_mach(self):
        inst = rfr(self.gamma, 1.5)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_t_t_star(self):
        inst = rfr(self.gamma, t_t_st=0.7525039918710987, flow_type=FS.SUPER_SONIC)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_p_p_star(self):
        inst = rfr(self.gamma, p_p_st=0.5783132530120482)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_rho_rho_star(self):
        inst = rfr(self.gamma, rho_rho_st=0.7685185185185186)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_po_po_star(self):
        inst = rfr(self.gamma, po_po_st=1.1215452274944158, flow_type=FS.SUPER_SONIC)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_to_to_star(self):
        inst = rfr(self.gamma, to_to_st=0.9092756568442442, flow_type=FS.SUPER_SONIC)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_from_u_u_star(self):
        inst = rfr(self.gamma, u_u_st=1.3012048192771082)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.4869, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(0.7593, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(0.586, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(0.7718, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.1152, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.2957, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.choked_flow

    def test_rayleigh_choked_flow(self):
        inst = rfr(self.gamma, u_u_st=1.3012048192771082)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flow_type == FS.SUPER_SONIC
        assert inst.t_t_st == approx(0.7525, rel=1e-4)
        assert inst.p_p_st == approx(0.5783, rel=1e-4)
        assert inst.rho_rho_st == approx(0.7685, rel=1e-2)
        assert inst.po_po_st == approx(1.1215, rel=1e-4)
        assert inst.to_to_st == approx(0.9093, rel=1e-4)
        assert inst.u_u_st == approx(1.3012, rel=1e-4)
        assert inst.choked_heat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gas_constant_r == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwn_strm_mach == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_t_t_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_p_p_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_rho_rho_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_po_po_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_u_u_st == approx(1.0, rel=1e-4)
        assert inst.dwn_strm_to_to_st == approx(1.0, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwn_strm_t_t_st / inst.t_t_st, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwn_strm_p_p_st / inst.p_p_st, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwn_strm_rho_rho_st / inst.rho_rho_st, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwn_strm_po_po_st / inst.po_po_st, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwn_strm_to_to_st / inst.to_to_st, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwn_strm_u_u_st / inst.u_u_st, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert inst.choked_flow
