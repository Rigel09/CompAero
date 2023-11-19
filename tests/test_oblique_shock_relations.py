from math import radians

import pytest
from pytest import approx

from CompAero.internal import InvalidOptionCombinationError, ShockType
from CompAero.oblique_shock_relations import ObliqueShockRelations as osr


class TestObliqueShockClassFuncs:
    gamma = 1.4

    #######################################################################################
    # Test the Functions for Subsonic Case
    #######################################################################################

    # TODO: Figure out how to handle some of the functions that calculate mach from a value?
    # What to do for the subsonic case, raise an exception?
    # Can it be detected? I.E. value always greater than X for mach > 1?
    def test_subsonic_p2_p1(self):
        with pytest.raises(ValueError):
            osr.calc_mach_normal_ahead_shock(0.5, self.gamma)

    #######################################################################################
    # Test the Functions for Superonic Case
    #######################################################################################
    # TODO: This function needs more tests and more checking on output validation
    def test_supersonic_calc_mach_normal_1(self):
        assert osr.calc_mach_normal_ahead_shock(1.5, radians(53.61)) == approx(1.2075, rel=1e-4)

    def test_supersonic_calc_mach_from_mach_normal_1(self):
        assert osr.calc_mach_ahead_shock_from_mach_normal_ahead_shock(
            1.207496034, radians(53.61)
        ) == approx(1.5, rel=1e-1)

    def test_supersonic_calc_beta_from_mach_normal_1(self):
        assert osr.calc_beta_from_mach_mach_normal_ahead_shock(1.5, 1.25340) == approx(
            radians(56.67868), rel=1e-4
        )

    def test_supersonic_calc_mach_2(self):
        assert osr.calc_mach_behind_shock(0.83754, radians(8.56287), radians(53.61)) == approx(
            1.18349, rel=1e-4
        )

    def test_supersonic_calc_mn2_from_m2(self):
        assert osr.calc_mach_normal_behind_shock_from_mach_behind_shock(
            1.11438, radians(56.67868), radians(10.0)
        ) == approx(0.81073, rel=1e-4)

    def test_supersonic_calc_theta_from_beta_mach(self):
        assert osr.calc_theta_from_theta_beta_mach(radians(56.67868), 1.5, self.gamma) == approx(
            radians(10.0), rel=1e-1
        )

    def test_supersonic_calc_beta_from_theta_mach_weak(self):
        assert osr.calc_beta_from_theta_beta_mach_weak(radians(10.0), 1.5, self.gamma) == approx(
            radians(56.67868), rel=1e-4
        )

    def test_supersonic_calc_beta_from_theta_mach_strong(self):
        assert osr.calc_beta_from_theta_beta_mach_strong(radians(10.0), 1.5, self.gamma) == approx(
            radians(75.99487), rel=1e-4
        )

    def test_supersonic_calc_mach_from_theta_beta(self):
        assert osr.calc_mach_from_theta_beta_mach(
            radians(56.67868), radians(10.0), self.gamma
        ) == approx(1.5, rel=1e-4)

    def test_supersonic_max_flow_deflection_angle(self):
        assert osr.calc_max_flow_deflection_angle(radians(66.5888), 1.5, self.gamma) == approx(
            radians(12.11267), rel=1e-4
        )

    def test_supersonic_max_shock_angle(self):
        assert osr.calc_max_shock_angle(1.5, self.gamma) == approx(radians(66.5888), rel=1e-4)

    def test_supersonic_calc_mach_wave_angle(self):
        assert osr.calc_mach_wave_angle(1.5) == approx(radians(41.8103149), rel=1e-5)

    def test_supersonic_calc_mach_from_mach_wave_angle(self):
        assert osr.calc_mach_from_mach_wave_angle(radians(41.8103149)) == approx(1.5, rel=1e-1)

    # TODO: Figure out someway to test the plotting feature of the TBM chart


###########################################################################################################
# Test the class construction
###########################################################################################################
class TestObliqueShockRelationsClass:
    gamma = 1.4

    def test_construction_from_mach_flow_deflection(self):
        inst = osr(self.gamma, mach=1.5, wedge_angle=10.0)
        assert inst.shock_angle == approx(56.67868, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.25340, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.81073, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(1.66619, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.43450, rel=1e-4)
        assert inst.t2_t1 == approx(1.16151, rel=1e-4)
        assert inst.po2_po1 == approx(0.98660, rel=1e-4)
        assert inst.po2_p1 == approx(2.56720)
        assert inst.mach2 == approx(1.11438, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_mach_norm1_shock_angle(self):
        inst = osr(self.gamma, mn1=1.2534, shock_angle=56.67868)
        assert inst.shock_angle == approx(56.67868, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.25340, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.81073, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(1.66619, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.43450, rel=1e-4)
        assert inst.t2_t1 == approx(1.16151, rel=1e-4)
        assert inst.po2_po1 == approx(0.98660, rel=1e-4)
        assert inst.po2_p1 == approx(2.5672, rel=1e-4)
        assert inst.mach2 == approx(1.11438, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_shock_angle_wedge_angle(self):
        inst = osr(self.gamma, wedge_angle=10.0, shock_angle=56.67868)
        assert inst.shock_angle == approx(56.67868, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.25340, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.81073, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(1.66619, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.43450, rel=1e-4)
        assert inst.t2_t1 == approx(1.16151, rel=1e-4)
        assert inst.po2_po1 == approx(0.98660, rel=1e-4)
        assert inst.po2_p1 == approx(2.5672, rel=1e-4)
        assert inst.mach2 == approx(1.11438, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_shock_angle_wedge_angle_m2_fail(self):
        with pytest.raises(ValueError):
            inst = osr(self.gamma, wedge_angle=10.0, shock_angle=56.67868, m2=2.5672)

    def test_construction_from_shock_angle_wedge_angle_m2(self):
        inst = osr(self.gamma, wedge_angle=10.0, shock_angle=23.01624, m2=3.13545)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_shock_angle_wedge_angle(self):
        inst = osr(self.gamma, wedge_angle=10.0, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_po2_p1(self):
        inst = osr(self.gamma, po2_p1=3.35977, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_p2_p1(self):
        inst = osr(self.gamma, p2_p1=2.40876, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_rho2_rho1(self):
        inst = osr(self.gamma, rho2_rho1=1.83768, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_t2_t1(self):
        inst = osr(self.gamma, t2_t1=1.31077, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_po2_po1(self):
        inst = osr(self.gamma, po2_po1=0.93423, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_from_mn2_shock_angle(self):
        inst = osr(self.gamma, mn2=0.70619, shock_angle=23.01624)
        assert inst.shock_angle == approx(23.01624, rel=1e-4)
        assert inst.wedge_angle == approx(10.0, rel=1e-4)
        assert inst.mach_normal_1 == approx(1.48577, rel=1e-4)
        assert inst.mach_normal_2 == approx(0.70619, rel=1e-4)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(3.8, rel=1e-1)
        assert inst.p2_p1 == approx(2.40876, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.83768, rel=1e-4)
        assert inst.t2_t1 == approx(1.31077, rel=1e-4)
        assert inst.po2_po1 == approx(0.93423, rel=1e-4)
        assert inst.po2_p1 == approx(3.35977, rel=1e-4)
        assert inst.mach2 == approx(3.13545, rel=1e-4)
        assert inst.shock_type == ShockType.WEAK

    def test_construction_not_enough_args(self):
        with pytest.raises(InvalidOptionCombinationError):
            osr(self.gamma)
