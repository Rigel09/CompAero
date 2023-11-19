from pytest import approx, raises

from CompAero.rocket_nozzle import (
    max_thrust_coefficient,
    min_thrust_coefficient,
    thrust_coefficient,
)


class Test_Rocket_Nozzle:
    gamma = 1.25

    def test_thrust_coefficient(self) -> None:
        assert thrust_coefficient(self.gamma, 4.0, 0.03952, 0.001) == approx(1.589749, rel=1e-4)
        assert thrust_coefficient(self.gamma, 5.0, 0.03952, 0.001) == approx(1.628265, rel=1e-4)
        assert thrust_coefficient(self.gamma, 10.0, 0.03952, 0.001) == approx(1.820845, rel=1e-4)
        assert thrust_coefficient(self.gamma, 1.0, 0.03952, 0.001) == approx(1.474202, rel=1e-4)

    def test_invalid_thrust_coefficient(self) -> None:
        with raises(ValueError):
            thrust_coefficient(1.4, 0.9, 0.003, 0.001)

    def test_max_thrust_coefficient(self) -> None:
        assert max_thrust_coefficient(self.gamma, 4.0, 0.03952) == approx(1.435686, rel=1e-4)
        assert max_thrust_coefficient(1.5, 2.0, 0.03952) == approx(1.394496, rel=1e-4)
        assert max_thrust_coefficient(1.7, 4.0, 0.03952) == approx(1.38162, rel=1e-4)
        assert max_thrust_coefficient(1.9, 4.0, 0.03952) == approx(1.37786, rel=1e-4)

    def test_invalid_max_thrust_coefficient(self) -> None:
        with raises(ValueError):
            max_thrust_coefficient(1.4, 0.9, 0.003)

    def test_min_thrust_coefficient(self) -> None:
        assert min_thrust_coefficient(1.0) == approx(0.1843, rel=1e-4)
        assert min_thrust_coefficient(4.0) == approx(0.836842, rel=1e-4)
        assert min_thrust_coefficient(6.0) == approx(0.99536, rel=1e-4)
        assert min_thrust_coefficient(10.0) == approx(1.174262, rel=1e-4)

    def test_invalid_min_thrust_coefficient(self) -> None:
        with raises(ValueError):
            min_thrust_coefficient(0.9)
