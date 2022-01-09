from pytest import approx
from CompAero.RayleighFlowRelations import RayleighFlowRelations as rfr
from CompAero.internal import FlowState as FS


class TestRayleighClassFuncs:
    gamma = 1.4

    # Test the Functions for Subsonic Case
    #######################################################################################
    def test_subsonic_t_tstar(self):
        assert rfr.calc_T_Tstar(0.5, self.gamma) == approx(0.7901, rel=1e-4)

    def test_subsonic_mach_from_t_tstar(self):
        assert rfr.calc_mach_from_T_TStar(0.79012345, self.gamma, FS.SUB_SONIC) == approx(0.5, rel=1e-2)

    def test_subsonic_p_pstar(self):
        assert rfr.calc_P_Pstar(0.5, self.gamma) == approx(1.7778, rel=1e-4)

    def test_subsonic_mach_from_p_pstar(self):
        assert rfr.calc_mach_from_P_PStar(1.77777777, self.gamma) == approx(0.5, rel=1e-2)

    def test_subsonic_rho_rhoStar(self):
        assert rfr.calc_Rho_RhoStar(0.5, self.gamma) == approx(2.25, rel=1e-4)

    def test_subsonic_mach_from_rho_rhoStar(self):
        assert rfr.calc_mach_from_Rho_RhoStar(2.25, 1.4) == approx(0.5, 1e-3)

    def test_subsonic_p0_p0Star(self):
        assert rfr.calc_Po_PoStar(0.5, self.gamma) == approx(1.1141, rel=1e-4)

    def test_subsonic_mach_from_p0_p0Star(self):
        assert rfr.calc_mach_from_Po_PoStar(1.11405250, self.gamma, flowType=FS.SUB_SONIC) == approx(
            0.5, 1e-3
        )

    def test_subsonic_t0_t0Star(self):
        assert rfr.calc_To_ToSt(0.5, self.gamma) == approx(0.6914, rel=1e-4)

    def test_subsonic_mach_from_t0_t0Star(self):
        assert rfr.calc_mach_from_To_ToStar(0.69135802, self.gamma, flowType=FS.SUB_SONIC) == approx(
            0.5, rel=1e-3
        )

    def test_subsonic_u_uStar(self):
        assert rfr.calc_U_UStar(0.5, self.gamma) == approx(0.44444, rel=1e-4)

    def test_subsonic_mach_from_u_uStar(self):
        assert rfr.calc_mach_from_U_USt(0.44444444, self.gamma) == approx(0.5, rel=1e-3)

    # Test the Functions for Supersonic Case
    #######################################################################################
    def test_supersonic_t_tstar(self):
        assert rfr.calc_T_Tstar(1.5, self.gamma) == approx(0.7525, rel=1e-4)

    def test_supersonic_mach_from_t_tstar(self):
        assert rfr.calc_mach_from_T_TStar(0.75250399, self.gamma, flowType=FS.SUPER_SONIC) == approx(
            1.5, rel=1e-2
        )

    def test_supersonic_p_pstar(self):
        assert rfr.calc_P_Pstar(1.5, self.gamma) == approx(0.5783, rel=1e-4)

    def test_supersonic_mach_from_p_pstar(self):
        assert rfr.calc_mach_from_P_PStar(0.57831325, self.gamma) == approx(1.5, rel=1e-2)

    def test_supersonic_rho_rhoStar(self):
        assert rfr.calc_Rho_RhoStar(1.5, self.gamma) == approx(0.7685, rel=1e-4)

    def test_supersonic_mach_from_rho_rhoStar(self):
        assert rfr.calc_mach_from_Rho_RhoStar(0.7685185185185186, 1.4) == approx(1.5, 1e-3)

    def test_supersonic_p0_p0Star(self):
        assert rfr.calc_Po_PoStar(1.5, self.gamma) == approx(1.1215, rel=1e-4)

    def test_supersonic_mach_from_p0_p0Star(self):
        assert rfr.calc_mach_from_Po_PoStar(1.12154522, self.gamma) == approx(1.5, 1e-3)

    def test_supersonic_t0_t0Star(self):
        assert rfr.calc_To_ToSt(1.5, self.gamma) == approx(0.9093, rel=1e-4)

    def test_supersonic_mach_from_t0_t0Star(self):
        assert rfr.calc_mach_from_To_ToStar(0.90927565, self.gamma, flowType=FS.SUPER_SONIC) == approx(
            1.5, rel=1e-3
        )

    def test_supersonic_u_uStar(self):
        assert rfr.calc_U_UStar(1.5, self.gamma) == approx(1.3012, rel=1e-4)

    def test_supersonic_mach_from_u_uStar(self):
        assert rfr.calc_mach_from_U_USt(1.30120481, self.gamma) == approx(1.5, rel=1e-3)


class TestRayleighClassSubsonic:
    gamma = 1.4

    def test_rayleigh_from_mach(self):
        inst = rfr(self.gamma, 0.5)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_t_tStar(self):
        inst = rfr(self.gamma, t_tSt=0.7901234567901234, flowType=FS.SUB_SONIC)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_p_pStar(self):
        inst = rfr(self.gamma, p_pSt=1.7777777777777777)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_rho_rhoStar(self):
        inst = rfr(self.gamma, rho_rhoSt=2.250)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_po_poStar(self):
        inst = rfr(self.gamma, po_poSt=1.114052503180089, flowType=FS.SUB_SONIC)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_to_toStar(self):
        inst = rfr(self.gamma, to_toSt=0.691358024691358, flowType=FS.SUB_SONIC)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_u_uStar(self):
        inst = rfr(self.gamma, u_uSt=0.4444444444444444)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(0.7525, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0148, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.3388, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.3192, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0295, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(0.7580, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9415, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_choked_flow(self):
        inst = rfr(self.gamma, u_uSt=0.4444444444444444)
        inst.simulate_heat_addition(200000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert inst.flowType == FS.SUB_SONIC
        assert inst.t_tSt == approx(0.7901, rel=1e-4)
        assert inst.p_pSt == approx(1.7778, rel=1e-4)
        assert inst.rho_rhoSt == approx(2.25, rel=1e-2)
        assert inst.po_poSt == approx(1.1141, rel=1e-4)
        assert inst.to_toSt == approx(0.6914, rel=1e-4)
        assert inst.u_uSt == approx(0.44444, rel=1e-4)
        assert inst.chokedHeat == approx(123410.0, rel=1e-1)
        assert inst.heat == approx(200000)
        assert inst.gasConstantR == approx(275.2, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(1.0, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(474.304, rel=1e-3)
        assert inst.chockedFlow


class TestRayleighClassSuperSonic:
    gamma = 1.4

    def test_rayleigh_from_mach(self):
        inst = rfr(self.gamma, 1.5)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_t_tStar(self):
        inst = rfr(self.gamma, t_tSt=0.7525039918710987, flowType=FS.SUPER_SONIC)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_p_pStar(self):
        inst = rfr(self.gamma, p_pSt=0.5783132530120482)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_rho_rhoStar(self):
        inst = rfr(self.gamma, rho_rhoSt=0.7685185185185186)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_po_poStar(self):
        inst = rfr(self.gamma, po_poSt=1.1215452274944158, flowType=FS.SUPER_SONIC)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_to_toStar(self):
        inst = rfr(self.gamma, to_toSt=0.9092756568442442, flowType=FS.SUPER_SONIC)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_from_u_uStar(self):
        inst = rfr(self.gamma, u_uSt=1.3012048192771082)
        inst.simulate_heat_addition(1000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(1000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.4869, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(0.7593, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(0.586, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(0.7718, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.1152, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.2957, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(0.9126, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(276.1955, rel=1e-3)
        assert not inst.chockedFlow

    def test_rayleigh_choked_flow(self):
        inst = rfr(self.gamma, u_uSt=1.3012048192771082)
        inst.simulate_heat_addition(100000, 275.2, 287)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.flowType == FS.SUPER_SONIC
        assert inst.t_tSt == approx(0.7525, rel=1e-4)
        assert inst.p_pSt == approx(0.5783, rel=1e-4)
        assert inst.rho_rhoSt == approx(0.7685, rel=1e-2)
        assert inst.po_poSt == approx(1.1215, rel=1e-4)
        assert inst.to_toSt == approx(0.9093, rel=1e-4)
        assert inst.u_uSt == approx(1.3012, rel=1e-4)
        assert inst.chokedHeat == approx(27582.056, rel=1e-1)
        assert inst.heat == approx(100000)
        assert inst.gasConstantR == approx(287.0, rel=1e-1)
        assert inst.cp == approx(1004.5, rel=1e-1)
        assert inst.dwnStrmMach == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_t_tSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_p_pSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_rho_rhoSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_po_poSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_u_uSt == approx(1.0, rel=1e-4)
        assert inst.dwnStrm_to_toSt == approx(1.0, rel=1e-4)
        assert inst.t2_t1 == approx(inst.dwnStrm_t_tSt / inst.t_tSt, rel=1e-4)
        assert inst.p2_p1 == approx(inst.dwnStrm_p_pSt / inst.p_pSt, rel=1e-4)
        assert inst.rho2_rho1 == approx(inst.dwnStrm_rho_rhoSt / inst.rho_rhoSt, rel=1e-4)
        assert inst.po2_po1 == approx(inst.dwnStrm_po_poSt / inst.po_poSt, rel=1e-4)
        assert inst.to2_to1 == approx(inst.dwnStrm_to_toSt / inst.to_toSt, rel=1e-4)
        assert inst.u2_u1 == approx(inst.dwnStrm_u_uSt / inst.u_uSt, rel=1e-4)
        assert inst.to1 == approx(275.2, rel=1e-1)
        assert inst.to2 == approx(374.752, rel=1e-3)
        assert inst.chockedFlow

