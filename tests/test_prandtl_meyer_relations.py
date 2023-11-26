# pylint: disable=missing-module-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring


import math

import pytest
from pytest import approx

from CompAero.internal import InvalidOptionCombinationError
from CompAero.prandtl_meyer import PrandtlMeyer as pm


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
    # pylint: disable=too-few-public-methods
    gamma = 1.4

    def test_subsonic_construction_given_mach(self):
        inst = pm(self.gamma, mach=0.5)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(0.5, rel=1e-1)
        assert math.isnan(inst.nu)
        assert math.isnan(inst.mu)
        assert math.isnan(inst.deflection_angle)
        assert math.isnan(inst.down_stream_nu)
        assert math.isnan(inst.down_stream_mu)
        assert math.isnan(inst.down_stream_mach)


class TestPrandtlMeyerClassSupersonic:
    gamma = 1.4

    def test_supersonic_construction_given_mach(self):
        inst = pm(self.gamma, mach=1.5)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert math.isnan(inst.deflection_angle)
        assert math.isnan(inst.down_stream_nu)
        assert math.isnan(inst.down_stream_mu)
        assert math.isnan(inst.down_stream_mach)

    def test_supersonic_construction_given_nu(self):
        inst = pm(self.gamma, nu=11.9052)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert math.isnan(inst.deflection_angle)
        assert math.isnan(inst.down_stream_nu)
        assert math.isnan(inst.down_stream_mu)
        assert math.isnan(inst.down_stream_mach)

    def test_supersonic_construction_given_mu(self):
        inst = pm(self.gamma, mu=41.81031)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert math.isnan(inst.deflection_angle)
        assert math.isnan(inst.down_stream_nu)
        assert math.isnan(inst.down_stream_mu)
        assert math.isnan(inst.down_stream_mach)

    def test_supersonic_construction_given_deflection_dwnstrm_mach(self):
        inst = pm(self.gamma, deflection_angle=10, down_stream_mach=1.84099)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflection_angle == approx(10)
        assert inst.down_stream_nu == approx(21.90521, rel=1e-4)
        assert inst.down_stream_mu == approx(32.9008, rel=1e-4)
        assert inst.down_stream_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_construction_given_deflection_dwnstrm_mu(self):
        inst = pm(self.gamma, deflection_angle=10, down_stream_mu=32.9008)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflection_angle == approx(10)
        assert inst.down_stream_nu == approx(21.90521, rel=1e-4)
        assert inst.down_stream_mu == approx(32.9008, rel=1e-4)
        assert inst.down_stream_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_construction_given_deflection_dwnstrm_nu(self):
        inst = pm(self.gamma, deflection_angle=10, down_stream_nu=21.90521)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflection_angle == approx(10)
        assert inst.down_stream_nu == approx(21.90521, rel=1e-4)
        assert inst.down_stream_mu == approx(32.9008, rel=1e-4)
        assert inst.down_stream_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_construction_given_defelction_angle_radians(self):
        inst = pm(
            self.gamma, deflection_angle=math.radians(10), in_degrees=False, down_stream_nu=21.90521
        )
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.nu == approx(11.9052, rel=1e-4)
        assert inst.mu == approx(41.81031, rel=1e-4)
        assert inst.deflection_angle == approx(10)
        assert inst.down_stream_nu == approx(21.90521, rel=1e-4)
        assert inst.down_stream_mu == approx(32.9008, rel=1e-4)
        assert inst.down_stream_mach == approx(1.84099, rel=1e-4)

    def test_supersonic_invalid_construction(self):
        with pytest.raises(InvalidOptionCombinationError):
            pm(self.gamma, deflection_angle=10)
