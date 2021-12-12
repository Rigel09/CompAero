from pytest import approx
from CompAero.IsentropecRelations import IsentropicRelations as isr


class TestIsentropicClassFuncs:
    gamma = 1.4

    # Test the Functions for Subsonic Case
    #######################################################################################
    def test_subsonic_func_to_t(self):
        assert isr.calc_T0_T(0.5, self.gamma) == approx(1.0500, rel=1e-4)

    def test_subsonic_func_mach_from_t0_t(self):
        assert isr.calcMachFrom_T0_T(1.05, self.gamma) == approx(0.5, rel=1e-3)

    def test_subsonic_func_p0_p(self):
        assert isr.calc_p0_p(0.5, self.gamma) == approx(1.1862, rel=1e-4)

    def test_subsonic_func_mach_from_p0_p(self):
        assert isr.calcMachFrom_p0_p(1.186212645674475, self.gamma) == approx(0.5, rel=1e-3)

    def test_subsonic_func_rho0_rho(self):
        assert isr.calc_rho0_rho(0.5, self.gamma) == approx(1.1297, rel=1e-4)

    def test_subsonic_func_mach_from_rho0_rho(self):
        assert isr.calcMachFrom_rho0_rho(1.1297263272993634, self.gamma) == approx(0.5, rel=1e-3)

    def test_subsonic_func_A_Astar(self):
        assert isr.calc_A_Astar(0.5, self.gamma) == approx(1.3398, rel=1e-4)

    def test_subsonic_func_mach_from_A_Astar(self):
        assert isr.calcMachFrom_A_Astar(1.33984375, self.gamma, flowType="Subsonic") == approx(0.5, rel=1e-4)

    # Test the Functions for Supersonic Case
    #######################################################################################
    def test_supersonic_func_to_t(self):
        assert isr.calc_T0_T(1.5, self.gamma) == approx(1.45, rel=1e-4)

    def test_supersonic_func_mach_from_t0_t(self):
        assert isr.calcMachFrom_T0_T(1.45, self.gamma) == approx(1.5, rel=1e-3)

    def test_supersonic_func_p0_p(self):
        assert isr.calc_p0_p(1.5, self.gamma) == approx(3.6710, rel=1e-4)

    def test_supersonic_func_mach_from_p0_p(self):
        assert isr.calcMachFrom_p0_p(3.671030714559521, self.gamma) == approx(1.5, rel=1e-3)

    def test_supersonic_func_rho0_rho(self):
        assert isr.calc_rho0_rho(1.5, self.gamma) == approx(2.5317, rel=1e-4)

    def test_supersonic_func_mach_from_rho0_rho(self):
        assert isr.calcMachFrom_rho0_rho(2.5317453011566733, self.gamma) == approx(1.5, rel=1e-3)

    def test_supersonic_func_A_Astar(self):
        assert isr.calc_A_Astar(1.5, self.gamma) == approx(1.1762, rel=1e-4)

    def test_supersonic_func_mach_from_A_Astar(self):
        assert isr.calcMachFrom_A_Astar(1.17616705, self.gamma, flowType="Supersonic") == approx(
            1.5, rel=1e-4
        )


class TestIsentropicClassSubsonic:
    gamma = 1.4

    def test_subsonic_construction_from_mach(self):
        inst = isr(gamma=self.gamma, mach=0.5)
        assert inst.mach == approx(0.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(1.1862, rel=1e-4)
        assert inst.t0_t == approx(1.05, rel=1e-3)
        assert inst.rho0_rho == approx(1.1297, rel=1e-4)
        assert inst.a_aStar == approx(1.3398, rel=1e-4)
        assert inst.flowType == "Subsonic"

    def test_subsonic_construction_from_p0_p(self):
        inst = isr(gamma=self.gamma, p0_p=1.186212645674475)
        assert inst.mach == approx(0.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(1.1862, rel=1e-4)
        assert inst.t0_t == approx(1.05, rel=1e-3)
        assert inst.rho0_rho == approx(1.1297, rel=1e-4)
        assert inst.a_aStar == approx(1.3398, rel=1e-4)
        assert inst.flowType == "Subsonic"

    def test_subsonic_construction_from_t0_t(self):
        inst = isr(gamma=self.gamma, t0_t=1.05)
        assert inst.mach == approx(0.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(1.1862, rel=1e-4)
        assert inst.t0_t == approx(1.05, rel=1e-3)
        assert inst.rho0_rho == approx(1.1297, rel=1e-4)
        assert inst.a_aStar == approx(1.3398, rel=1e-4)
        assert inst.flowType == "Subsonic"

    def test_subsonic_construction_from_rho0_rho(self):
        inst = isr(gamma=self.gamma, rho0_rho=1.1297263272993634)
        assert inst.mach == approx(0.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(1.1862, rel=1e-4)
        assert inst.t0_t == approx(1.05, rel=1e-3)
        assert inst.rho0_rho == approx(1.1297, rel=1e-4)
        assert inst.a_aStar == approx(1.3398, rel=1e-4)
        assert inst.flowType == "Subsonic"

    def test_subsonic_construction_from_A_Astar(self):
        inst = isr(gamma=self.gamma, a_aStar=1.33984375, flowType="Subsonic")
        assert inst.mach == approx(0.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(1.1862, rel=1e-4)
        assert inst.t0_t == approx(1.05, rel=1e-3)
        assert inst.rho0_rho == approx(1.1297, rel=1e-4)
        assert inst.a_aStar == approx(1.3398, rel=1e-4)
        assert inst.flowType == "Subsonic"


class TestIsentropicClassSupersonic:
    gamma = 1.4

    def test_supersonic_construction_from_mach(self):
        inst = isr(gamma=self.gamma, mach=1.5)
        assert inst.mach == approx(1.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(3.6710, rel=1e-4)
        assert inst.t0_t == approx(1.45, rel=1e-3)
        assert inst.rho0_rho == approx(2.5317, rel=1e-4)
        assert inst.a_aStar == approx(1.1762, rel=1e-4)
        assert inst.flowType == "Supersonic"

    def test_supersonic_construction_from_p0_p(self):
        inst = isr(gamma=self.gamma, p0_p=3.671030714559521)
        assert inst.mach == approx(1.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(3.6710, rel=1e-4)
        assert inst.t0_t == approx(1.45, rel=1e-3)
        assert inst.rho0_rho == approx(2.5317, rel=1e-4)
        assert inst.a_aStar == approx(1.1762, rel=1e-4)
        assert inst.flowType == "Supersonic"

    def test_supersonic_construction_from_t0_t(self):
        inst = isr(gamma=self.gamma, t0_t=1.45)
        assert inst.mach == approx(1.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(3.6710, rel=1e-4)
        assert inst.t0_t == approx(1.45, rel=1e-3)
        assert inst.rho0_rho == approx(2.5317, rel=1e-4)
        assert inst.a_aStar == approx(1.1762, rel=1e-4)
        assert inst.flowType == "Supersonic"

    def test_supersonic_construction_from_rho0_rho(self):
        inst = isr(gamma=self.gamma, rho0_rho=2.5317453011566733)
        assert inst.mach == approx(1.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(3.6710, rel=1e-4)
        assert inst.t0_t == approx(1.45, rel=1e-3)
        assert inst.rho0_rho == approx(2.5317, rel=1e-4)
        assert inst.a_aStar == approx(1.1762, rel=1e-4)
        assert inst.flowType == "Supersonic"

    def test_supersonic_construction_from_A_Astar(self):
        inst = isr(gamma=self.gamma, a_aStar=1.17616705)
        assert inst.mach == approx(1.5, rel=1e-2)
        assert inst.gamma == approx(self.gamma, rel=1e-2)
        assert inst.p0_p == approx(3.6710, rel=1e-4)
        assert inst.t0_t == approx(1.45, rel=1e-3)
        assert inst.rho0_rho == approx(2.5317, rel=1e-4)
        assert inst.a_aStar == approx(1.1762, rel=1e-4)
        assert inst.flowType == "Supersonic"

