from numpy import pi, sin, cos
import pytest
import quadpy


@pytest.mark.parametrize('k', range(1, 6))
def test_simple(k):
    val, _ = quadpy.triangle.adaptive_integrate(
            lambda x: sin(k * pi * x[0]) * sin(k * pi * x[1]),
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            1.0e-10
            )
    exact = (2 - pi*k*sin(pi*k) - 2*cos(pi*k)) / (2 * pi**2 * k**2)

    assert abs(exact - val) < 1.0e-10
    return


def test_vector(k=1):
    val, _ = quadpy.triangle.adaptive_integrate(
            lambda x: [
                sin(k * pi * x[0]) * sin(k * pi * x[1]),
                cos(k * pi * x[0]) * cos(k * pi * x[1]),
                ],
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            1.0e-10
            )
    exact = [
        (2 - pi*k*sin(pi*k) - 2*cos(pi*k)) / (2 * pi**2 * k**2),
        sin(pi * k) / (2*pi*k)
        ]

    assert (abs(exact - val) < 1.0e-10).all()
    return


if __name__ == '__main__':
    test_vector(1.1)
