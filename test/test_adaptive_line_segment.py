from numpy import pi, sin, cos, array
import pytest
import quadpy


def test_simple():
    val, _ = quadpy.line_segment.integrate_adaptive(sin, [0.0, pi], 1.0e-10)
    exact = 2.0
    assert abs(exact - val) < 1.0e-10

    val, _ = quadpy.line_segment.integrate_adaptive(
        lambda x: x * sin(x), [0.0, pi], 1.0e-10
    )
    exact = pi
    assert abs(exact - val) < 1.0e-10


@pytest.mark.parametrize("k", range(1, 6))
def test_vector_valued(k):
    val, _ = quadpy.line_segment.integrate_adaptive(
        lambda x: array([x * sin(k * x), x * cos(k * x)]), [0.0, pi], 1.0e-10
    )
    exact = [
        (sin(pi * k) - pi * k * cos(pi * k)) / k ** 2,
        (cos(pi * k) + pi * k * sin(pi * k) - 1.0) / k ** 2,
    ]
    assert (abs(exact - val) < 1.0e-10).all()


def test_predefined_intervals():
    k = 5
    val, _ = quadpy.line_segment.integrate_adaptive(
        lambda x: x * sin(k * x),
        [[0.0, 0.3 * pi, 0.5 * pi], [0.3 * pi, 0.5 * pi, pi]],
        1.0e-10,
    )
    exact = (sin(pi * k) - pi * k * cos(pi * k)) / k ** 2
    assert abs(exact - val) < 1.0e-10


@pytest.mark.parametrize("k", range(4, 12))
def test_sink(k):
    val, _ = quadpy.line_segment.integrate_adaptive(
        lambda x: sin(k * x), [0.0, pi], 1.0e-10
    )
    exact = (1.0 - cos(k * pi)) / k
    assert abs(exact - val) < 1.0e-10
    return


if __name__ == "__main__":
    # test_predefined_intervals()
    test_vector_valued(1)
    # test_simple()
    # test_sink(5)
