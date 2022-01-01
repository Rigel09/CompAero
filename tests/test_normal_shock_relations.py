from math import gamma
from pytest import approx
import pytest
from CompAero.NormalShockRelations import NormalShockRelations as nsr


class TestNormalShockClassFuncs:
    gamma = 1.4

    #######################################################################################
    # Test the Functions for Subsonic Case
    #######################################################################################

    # TODO: Figure out how to handle some of the functions that calculate mach from a value?
    # What to do for the subsonic case, raise an exception?
    # Can it be detected? I.E. value always greater than X for mach > 1?

    def test_subsonic_p2_p1(self):
        with pytest.raises(ValueError):
            nsr.calc_p2_p1(0.5, self.gamma)

    def test_subsonic_rho2_rho1(self):
        with pytest.raises(ValueError):
            nsr.calc_rho2_rho1(0.5, self.gamma)

    def test_subsonic_t2_t1(self):
        with pytest.raises(ValueError):
            nsr.calc_T2_T1(0.5, self.gamma)

    def test_subsonic_downstream_mach(self):
        with pytest.raises(ValueError):
            nsr.calc__mach2(0.5, self.gamma)

    def test_subsonic_po2_po1(self):
        with pytest.raises(ValueError):
            nsr.calc_po2_po1(0.5, self.gamma)

    def test_subsonic_po2_p1(self):
        with pytest.raises(ValueError):
            nsr.calc_p2_p1(0.5, self.gamma)

    #######################################################################################
    # Test the Functions for Supersonic Case
    #######################################################################################
    def test_supersonic_p2_p1(self):
        assert nsr.calc_p2_p1(1.5, self.gamma) == approx(2.45833, rel=1e-4)

    def test_supersonic_mach_from_p2_p1(self):
        assert nsr.calcMachFrom_p2_p1(2.45833, self.gamma) == approx(1.5, rel=1e-1)

    def test_supersonic_rho2_rho1(self):
        assert nsr.calc_rho2_rho1(1.5, self.gamma) == approx(1.86207, rel=1e-4)

    def test_supersonic_mach_from_rho2_rho1(self):
        assert nsr.calcMachFrom_rho2_rho1(1.86207, self.gamma) == approx(1.5, rel=1e-1)

    def test_supersonic_t2_t1(self):
        assert nsr.calc_T2_T1(1.5, self.gamma) == approx(1.32022, rel=1e-4)

    def test_supersonic_mach_from_t2_t1(self):
        assert nsr.calcMachFrom_T2_T1(1.32022, self.gamma) == approx(1.5, rel=1e-1)

    def test_supersonic_downstream_mach(self):
        assert nsr.calc__mach2(1.5, self.gamma) == approx(0.70109, rel=1e-4)

    def test_supersonic_mach_from_downstream_mach(self):
        assert nsr.calcMachFrom_mach2(0.70109, self.gamma) == approx(1.5, rel=1e-4)

    def test_supersonic_po2_po1(self):
        assert nsr.calc_po2_po1(1.5, self.gamma) == approx(0.92979, rel=1e-4)

    def test_supersonic_mach_from_po2_po1(self):
        assert nsr.calcMachFrom_po2_po1(0.92979, self.gamma) == approx(1.5, rel=1e-4)

    def test_supersonic_po2_p1(self):
        assert nsr.calc_po2_p1(1.5, self.gamma) == approx(3.41327, rel=1e-4)

    def test_supersonic_mach_from_po2_p1(self):
        assert nsr.calcMachFrom_po2_p1(3.41327, self.gamma) == approx(1.5, rel=1e-1)


#######################################################################################
# Test the Different Construction of a class
#######################################################################################
class TestNormalShockRelationsConstructions:
    gamma = 1.4

    def test_construction_from_mach(self):
        inst = nsr(self.gamma, mach=1.5)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

    def test_construction_from_p2_p1(self):
        inst = nsr(self.gamma, p2_p1=2.45833)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

    def test_construction_from_rho2_rho1(self):
        inst = nsr(self.gamma, rho2_rho1=1.86207)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

    def test_construction_from_t2_t1(self):
        inst = nsr(self.gamma, t2_t1=1.32022)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

    def test_construction_from_po2_po1(self):
        inst = nsr(self.gamma, po2_po1=0.92979)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

    def test_construction_from_po2_p1(self):
        inst = nsr(self.gamma, po2_p1=3.41327)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

    def test_construction_from_mach2(self):
        inst = nsr(self.gamma, m2=0.70109)
        assert inst.gamma == approx(self.gamma, rel=1e-1)
        assert inst.mach == approx(1.5, rel=1e-1)
        assert inst.p2_p1 == approx(2.45833, rel=1e-4)
        assert inst.rho2_rho1 == approx(1.86207, rel=1e-4)
        assert inst.t2_t1 == approx(1.32022, rel=1e-4)
        assert inst.po2_po1 == approx(0.92979, rel=1e-4)
        assert inst.po2_p1 == approx(3.41327, rel=1e-4)
        assert inst.mach2 == approx(0.70109, rel=1e-4)

