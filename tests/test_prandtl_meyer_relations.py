import math
from numpy import isnat
from pytest import approx
import pytest
from CompAero.PrandtlMeyer import PrandtlMeyer as pm
from CompAero.internal import InvalidOptionCombinationError


#########################################################################################
# Test the static functions of the class
#########################################################################################
class TestPrandtlMeyerClassFuncs:
    gamma = 1.4

    def test_subsonic_nu_from_mach(self):
        assert pm.calc_nu(0.5, self.gamma) == 0

    def test_supersonic_nu_from_mach(self):
        assert pm.calc_nu(1.5, self.gamma) == approx(11.9052, rel=1e-4)

    def test_subsonic_mach_from_nu(self):
        assert pm.calc_mach_from_nu(0, self.gamma) == 1.0

    def test_supersonic_mach_from_nu(self):
        assert pm.calc_mach_from_nu(11.9052, self.gamma) == approx(1.5, rel=1e-1)


#########################################################################################
# Test the different construction methods of the class
#########################################################################################


class TestPrandtlMeyerClassSubsonic:
    gamma = 1.4

    def test_subsonic_construction_given_mach(self):
        inst = pm(self.gamma, mach=0.5)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert math.isnan(inst.nu)
        assert math.isnan(inst.mu)
        assert math.isnan(inst.deflectionAngle)
        assert math.isnan(inst.dwmStrm_nu)
        assert math.isnan(inst.dwmStrm_mu)
        assert math.isnan(inst.dwmStrm_mach)


class TestPrandtlMeyerClassSupersonic:
    gamma = 1.4

    def test_supersonic_construction_given_mach(self):
        inst = pm(self.gamma, mach=1.5)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert math.isnan(inst.deflectionAngle)
        assert math.isnan(inst.dwmStrm_nu)
        assert math.isnan(inst.dwmStrm_mu)
        assert math.isnan(inst.dwmStrm_mach)

    def test_supersonic_construction_given_nu(self):
        inst = pm(self.gamma, nu=11.9052)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert math.isnan(inst.deflectionAngle)
        assert math.isnan(inst.dwmStrm_nu)
        assert math.isnan(inst.dwmStrm_mu)
        assert math.isnan(inst.dwmStrm_mach)

    def test_supersonic_construction_given_mu(self):
        inst = pm(self.gamma, mu=41.81031)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert math.isnan(inst.deflectionAngle)
        assert math.isnan(inst.dwmStrm_nu)
        assert math.isnan(inst.dwmStrm_mu)
        assert math.isnan(inst.dwmStrm_mach)

    def test_supersonic_construction_given_deflection_dwnstrm_mach(self):
        inst = pm(self.gamma, deflectionAngle=10, dwnStreamMach=1.84099)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflectionAngle == approx(10)
        assert inst.dwmStrm_nu == approx(21.90521, rel=1e-4)
        assert inst.dwmStrm_mu == approx(32.9008, rel=1e-4)
        assert inst.dwmStrm_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_construction_given_deflection_dwnstrm_mu(self):
        inst = pm(self.gamma, deflectionAngle=10, dwnStreamMu=32.9008)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflectionAngle == approx(10)
        assert inst.dwmStrm_nu == approx(21.90521, rel=1e-4)
        assert inst.dwmStrm_mu == approx(32.9008, rel=1e-4)
        assert inst.dwmStrm_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_construction_given_deflection_dwnstrm_nu(self):
        inst = pm(self.gamma, deflectionAngle=10, dwnstreamNu=21.90521)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflectionAngle == approx(10)
        assert inst.dwmStrm_nu == approx(21.90521, rel=1e-4)
        assert inst.dwmStrm_mu == approx(32.9008, rel=1e-4)
        assert inst.dwmStrm_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_construction_given_defelction_angle_radians(self):
        inst = pm(self.gamma, deflectionAngle=math.radians(10), inDegrees=False, dwnstreamNu=21.90521)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflectionAngle == approx(10)
        assert inst.dwmStrm_nu == approx(21.90521, rel=1e-4)
        assert inst.dwmStrm_mu == approx(32.9008, rel=1e-4)
        assert inst.dwmStrm_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_invalid_construction(self):
        with pytest.raises(InvalidOptionCombinationError):
            pm(self.gamma, deflectionAngle=10)
