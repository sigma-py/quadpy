import numpy as np
import orthopy
import pytest
from mpmath import mp

import quadpy


@pytest.mark.parametrize(
    "scheme,tol", [(quadpy.e1r.gauss_laguerre(n), 1.0e-14) for n in range(1, 10)]
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    evaluator = orthopy.e1r.Eval(scheme.points.T, "normal")

    degree = None
    for k in range(scheme.degree + 2):
        approximate = scheme.integrate(lambda x: next(evaluator))
        exact = evaluator.int_p0 if k == 0 else 0.0
        err = np.abs(approximate - exact)
        if np.any(err > tol):
            degree = k - 1
            break

    max_err = np.max(err)
    assert degree >= scheme.degree, (
        f"{scheme.name} -- observed: {degree}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e1r.gauss_laguerre(1)])
def test_show(scheme):
    scheme.show()


def test_laguerre_mpmath():
    mp.dps = 51
    scheme = quadpy.e1r.gauss_laguerre(2, mode="mpmath")

    tol = 1.0e-50

    x1 = 2 - mp.sqrt(2)
    x2 = 2 + mp.sqrt(2)
    assert (abs(scheme.points_symbolic - [x1, x2]) < tol).all()

    w1 = (2 + mp.sqrt(2)) / 4
    w2 = (2 - mp.sqrt(2)) / 4
    assert (abs(scheme.weights_symbolic - [w1, w2]) < tol).all()


def test_laguerre_generalized_mpmath():
    mp.dps = 51
    scheme = quadpy.e1r.gauss_laguerre(2, alpha=1, mode="mpmath")

    tol = 1.0e-50

    x1 = 3 - mp.sqrt(3)
    x2 = 3 + mp.sqrt(3)
    assert (abs(scheme.points_symbolic - [x1, x2]) < tol).all()

    w1 = 2 / ((-1 + mp.sqrt(3)) ** 2 * (1 + 2 / (-1 + mp.sqrt(3)) ** 2))
    w2 = 2 / ((-1 - mp.sqrt(3)) ** 2 * (1 + 2 / (-1 - mp.sqrt(3)) ** 2))
    assert (abs(scheme.weights_symbolic - [w1, w2]) < tol).all()


if __name__ == "__main__":
    scheme_ = quadpy.e1r.gauss_laguerre(3)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
