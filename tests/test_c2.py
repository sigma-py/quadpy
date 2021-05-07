import numpy as np
import orthopy
import pytest
from helpers import find_best_scheme

import quadpy

schemes = (
    list(quadpy.c2.schemes.values())
    + [quadpy.c2.product(quadpy.c1.midpoint())]
    + [quadpy.c2.product(quadpy.c1.trapezoidal())]
    + [quadpy.c2.product(quadpy.c1.gauss_legendre(k)) for k in range(1, 5)]
    + [quadpy.c2.product(quadpy.c1.newton_cotes_closed(k)) for k in range(1, 5)]
    + [quadpy.c2.product(quadpy.c1.newton_cotes_open(k)) for k in range(1, 6)]
)


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    # instantiate
    try:
        scheme = scheme()
    except TypeError:
        pass

    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    print(scheme)

    quad = quadpy.c2.rectangle_points([-1.0, +1.0], [-1.0, +1.0])

    evaluator = orthopy.cn.Eval(scheme.points)

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), quad)
        exact = evaluator.int_p0 * 4 if k == 0 else 0.0
        err = np.abs(approximate - exact)
        max_err = max(max_err, np.max(err))
        if np.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), quad)
            exact = 2.0 if i == 0 else 0.0
            err = np.abs(approximate - exact)
            max_err = max(max_err, np.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
        )


@pytest.mark.parametrize("scheme", [quadpy.c2.product(quadpy.c1.gauss_legendre(5))])
def test_show(scheme):
    scheme.show()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.c2.schemes.values(),
            degree,
            lambda pts: np.all((pts >= -1) & (pts <= 1)),
            lambda keys: len(keys - {"d4_a0", "d4_aa", "d4_ab", "zero2"}) == 0,
        )
        if best is None:
            break

        b = quadpy.c2.get_good_scheme(degree)

        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 22


if __name__ == "__main__":
    test_get_good_scheme()
