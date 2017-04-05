from numpy import pi, sin
import quadpy


def test_simple():
    val, _ = quadpy.triangle.adaptive_integrate(
            lambda x: sin(pi * x[0]) * sin(pi * x[1]),
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            1.0e-10
            )
    exact = 2 / pi**2

    assert abs(exact - val) < 1.0e-10


if __name__ == '__main__':
    test_simple()
