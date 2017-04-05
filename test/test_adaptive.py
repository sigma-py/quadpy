from numpy import pi, sin, cos, array
import pytest
import quadpy


def test_simple():
    val, _ = quadpy.line_segment.adaptive_integrate(
            sin, [0.0, pi], 1.0e-10
            )
    exact = 2.0
    assert abs(exact - val) < 1.0e-10

    val, _ = quadpy.line_segment.adaptive_integrate(
            lambda x: x * sin(x),
            [0.0, pi],
            1.0e-10
            )
    exact = pi
    assert abs(exact - val) < 1.0e-10


def test_vector_valued():
    val, _ = quadpy.line_segment.adaptive_integrate(
            lambda x: array([x * sin(x), x * cos(x)]),
            [0.0, pi],
            1.0e-10
            )
    exact = [pi, -2.0]
    assert (abs(exact - val) < 1.0e-10).all()


def test_predefined_intervals():
    k = 5
    val, _ = quadpy.line_segment.adaptive_integrate(
            lambda x: x * sin(k * x),
            [
                [0.0, 0.3*pi, 0.5*pi],
                [0.3*pi, 0.5*pi, pi],
            ],
            1.0e-10
            )
    exact = (sin(pi * k) - pi*k * cos(pi*k)) / k**2
    assert abs(exact - val) < 1.0e-10


@pytest.mark.parametrize('k', range(4, 12))
def test_sink(k):
    val, _ = quadpy.line_segment.adaptive_integrate(
            lambda x: sin(k*x),
            [0.0, pi],
            1.0e-10
            )
    exact = (1.0 - cos(k*pi)) / k
    assert abs(exact - val) < 1.0e-10
    return


if __name__ == '__main__':
    # test_predefined_intervals()
    test_vector_valued()
    # test_simple()
    # test_sink(5)
