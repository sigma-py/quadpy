import numpy as np
import orthopy
import pytest
from mpmath import mp

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.e1r2.gauss_hermite(n) for n in range(5, 12)]
    + [quadpy.e1r2.genz_keister(n) for n in range(8)],
)
def test_scheme(scheme, tol=1.0e-14):
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    evaluator = orthopy.e1r2.Eval(scheme.points.T, "physicists", "normal")

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


@pytest.mark.parametrize("scheme", [quadpy.e1r2.gauss_hermite(2)])
def test_show(scheme):
    scheme.show()


def test_hermite_mpmath():
    mp.dps = 51
    scheme = quadpy.e1r2.gauss_hermite(4, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.sqrt((3 - mp.sqrt(6)) / 2)
    x2 = mp.sqrt((3 + mp.sqrt(6)) / 2)
    assert (abs(scheme.points_symbolic - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = mp.sqrt(mp.pi) / 4 / (3 - mp.sqrt(6))
    w2 = mp.sqrt(mp.pi) / 4 / (3 + mp.sqrt(6))

    assert (abs(scheme.weights_symbolic - [w2, w1, w1, w2]) < tol).all()


if __name__ == "__main__":
    scheme_ = quadpy.e1r2.gauss_hermite(10)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
