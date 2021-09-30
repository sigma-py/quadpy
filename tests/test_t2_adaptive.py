import pytest
from numpy import cos, exp, pi, sin

import quadpy


@pytest.mark.parametrize("k", range(1, 6))
def test_simple(k):
    val, _ = quadpy.t2.integrate_adaptive(
        lambda x: sin(k * pi * x[0]) * sin(k * pi * x[1]),
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        1.0e-10,
    )
    exact = (2 - pi * k * sin(pi * k) - 2 * cos(pi * k)) / (2 * pi ** 2 * k ** 2)

    assert abs(exact - val) < 1.0e-10


def test_vector(k=1):
    triangles = [
        [[0.0, 0.0], [1.0, 0.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [[0.0, 1.0], [0.0, 1.0]],
    ]
    val, _ = quadpy.t2.integrate_adaptive(
        lambda x: [
            sin(k * pi * x[0]) * sin(k * pi * x[1]),
            cos(k * pi * x[0]) * cos(k * pi * x[1]),
            exp(k * pi * x[0]),
        ],
        triangles,
        1.0e-10,
    )
    exact = [
        (cos(pi * k) - 1) * (cos(pi * k) - 1) / (pi ** 2 * k ** 2),
        sin(pi * k) * sin(pi * k) / (pi ** 2 * k ** 2),
        (exp(pi * k) - 1) / (pi * k),
    ]

    assert (abs(exact - val) < 1.0e-10).all()


if __name__ == "__main__":
    # test_simple(1.1)
    test_vector(1.1)
