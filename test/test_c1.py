import math

import numpy as np
import orthopy
import pytest
from mpmath import mp

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.c1.midpoint()]
    + [quadpy.c1.trapezoidal()]
    + [quadpy.c1.clenshaw_curtis(k) for k in range(2, 10)]
    + [quadpy.c1.gauss_legendre(k) for k in range(1, 6)]
    + [quadpy.c1.gauss_lobatto(k) for k in range(2, 7)]
    + [quadpy.c1.gauss_kronrod(k) for k in range(2, 7)]
    + [quadpy.c1.gauss_patterson(k) for k in range(9)]
    + [quadpy.c1.gauss_radau(k) for k in range(2, 10)]
    + [quadpy.c1.fejer_1(k) for k in range(1, 10)]
    + [quadpy.c1.fejer_2(k) for k in range(1, 10)]
    + [quadpy.c1.newton_cotes_closed(k) for k in range(1, 5)]
    + [quadpy.c1.newton_cotes_open(k) for k in range(1, 5)],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    # https://github.com/nschloe/quadpy/issues/227
    assert scheme.weights.ndim == 1

    # test scheme.__str__
    print(scheme)

    degree = 0
    while True:
        # Set bounds such that the values are between 0.5 and 1.5.
        exact_val = 1.0 / (degree + 1)
        interval = np.array(
            [
                [0.5 ** (1.0 / (degree + 1)), 0.0, 0.0],
                [1.5 ** (1.0 / (degree + 1)), 0.0, 0.0],
            ]
        )
        interval = np.array([[0.3], [0.5]])
        val = scheme.integrate(lambda x: x[0] ** degree, interval)
        # same test with line embedded in R^2
        interval = np.array(
            [[0.5 ** (1.0 / (degree + 1)), 0.0], [1.5 ** (1.0 / (degree + 1)), 0.0]]
        )
        val = scheme.integrate(lambda x: x[0] ** degree, interval)
        if abs(exact_val - val) > 1.0e-12 * abs(exact_val):
            break
        if degree >= scheme.degree:
            break
        degree += 1
    assert degree == scheme.degree


@pytest.mark.parametrize(
    "scheme", [quadpy.c1.chebyshev_gauss_1(k) for k in range(1, 10)]
)
def test_cheb1_scheme(scheme):
    evaluator = orthopy.c1.chebyshev1.Eval(scheme.points, "normal")

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), [-1, 1])
        exact = math.sqrt(math.pi) if k == 0 else 0.0
        err = np.abs(approximate - exact)
        if np.any(err > 1.0e-14):
            break
        k += 1

    max_err = np.max(err)
    assert k - 1 >= scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize(
    "scheme", [quadpy.c1.chebyshev_gauss_2(k) for k in range(1, 10)]
)
def test_cheb2_scheme(scheme):
    evaluator = orthopy.c1.chebyshev2.Eval(scheme.points, "normal")

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), [-1, 1])
        exact = math.sqrt(math.pi / 2) if k == 0 else 0.0
        err = np.abs(approximate - exact)
        if np.any(err > 1.0e-14):
            break
        k += 1

    max_err = np.max(err)
    assert k - 1 >= scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.c1.newton_cotes_closed(5)])
def test_show(scheme):
    scheme.show()


def test_integrate_split():
    x = np.linspace(0.15, 0.702, 101)
    intervals = np.array([x[:-1], x[1:]])
    scheme = quadpy.c1.trapezoidal()
    val = scheme.integrate(
        lambda r: 0.5108
        / r ** 2
        / np.sqrt(2 * 1.158 + 2 / r - 0.5108 ** 2 / (2 * r ** 2)),
        intervals,
    )
    val = np.sum(val)
    reference = 0.961715
    assert abs(val - reference) < 1.0e-3 * reference


def test_legendre_mpmath():
    mp.dps = 50
    scheme = quadpy.c1.gauss_legendre(4, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.sqrt(mp.mpf(3) / 7 - mp.mpf(2) / 7 * mp.sqrt(mp.mpf(6) / 5))
    x2 = mp.sqrt(mp.mpf(3) / 7 + mp.mpf(2) / 7 * mp.sqrt(mp.mpf(6) / 5))
    assert (abs(scheme.points_symbolic - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = (18 + mp.sqrt(30)) / 36
    w2 = (18 - mp.sqrt(30)) / 36
    assert (abs(scheme.weights_symbolic - [w2, w1, w1, w2]) < tol).all()


def test_chebyshev1_sympy():
    scheme = quadpy.c1.chebyshev_gauss_1(4, mode="sympy")
    scheme_numpy = quadpy.c1.chebyshev_gauss_1(4, mode="numpy")

    flt = np.vectorize(float)
    tol = 1.0e-15

    assert (abs(flt(scheme.points) - scheme_numpy.points) < tol).all()
    assert (abs(flt(scheme.weights) - scheme_numpy.weights) < tol).all()


def test_chebyshev2_sympy():
    scheme = quadpy.c1.chebyshev_gauss_2(4, mode="sympy")
    scheme_numpy = quadpy.c1.chebyshev_gauss_2(4, mode="numpy")

    flt = np.vectorize(float)
    tol = 1.0e-15

    assert (abs(flt(scheme.points) - scheme_numpy.points) < tol).all()
    assert (abs(flt(scheme.weights) - scheme_numpy.weights) < tol).all()


def test_chebyshev1_mpmath():
    mp.dps = 50
    scheme = quadpy.c1.chebyshev_gauss_1(4, mode="mpmath")
    tol = 1.0e-50

    x1 = mp.cos(3 * mp.pi / 8)
    x2 = mp.cos(1 * mp.pi / 8)
    assert (abs(scheme.points_symbolic - [+x2, +x1, -x1, -x2]) < tol).all()

    w = mp.pi / 4
    tol = 1.0e-49
    assert (abs(scheme.weights_symbolic - [w, w, w, w]) < tol).all()


def test_chebyshev2_mpmath():
    mp.dps = 51
    scheme = quadpy.c1.chebyshev_gauss_2(4, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.cos(2 * mp.pi / 5)
    x2 = mp.cos(1 * mp.pi / 5)
    assert (abs(scheme.points_symbolic - [+x2, +x1, -x1, -x2]) < tol).all()

    w1 = mp.pi / 5 * mp.sin(2 * mp.pi / 5) ** 2
    w2 = mp.pi / 5 * mp.sin(1 * mp.pi / 5) ** 2
    assert (abs(scheme.weights_symbolic - [w2, w1, w1, w2]) < tol).all()


def test_jacobi_mpmath():
    mp.dps = 51
    scheme = quadpy.c1.gauss_jacobi(4, 1, 1, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.sqrt((7 - 2 * mp.sqrt(7)) / 21)
    x2 = mp.sqrt((7 + 2 * mp.sqrt(7)) / 21)
    assert (abs(scheme.points_symbolic - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = (5 + mp.sqrt(7)) / 15
    w2 = (5 - mp.sqrt(7)) / 15
    assert (abs(scheme.weights_symbolic - [w2, w1, w1, w2]) < tol).all()


def test_multidim():
    scheme = quadpy.c1.gauss_legendre(5)

    # simple scalar integration
    val = scheme.integrate(np.sin, [0.0, 1.0])
    assert val.shape == ()

    # scalar integration on 3 subdomains
    val = scheme.integrate(np.sin, [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]])
    assert val.shape == (3,)

    # scalar integration in 3D
    val = scheme.integrate(
        lambda x: x[0] + np.sin(x[1]) + np.cos(x[2]),
        [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
    )
    assert val.shape == ()

    # vector-valued integration on 3 subdomains
    val = scheme.integrate(
        lambda x: [np.sin(x), np.cos(x)], [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]]
    )
    assert val.shape == (2, 3)

    # vector-valued integration in 3D
    val = scheme.integrate(
        lambda x: [x[0] + np.sin(x[1]), np.cos(x[0]) * x[2]],
        [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
    )
    assert val.shape == (2,)

    # another vector-valued integration in 3D
    # This is one case where the integration routine may not properly recognize the
    # dimensionality of the domain. Use the `dim` parameter.
    val = scheme.integrate(
        lambda x: [
            x[0] + np.sin(x[1]),
            np.cos(x[0]) * x[2],
            np.sin(x[0]) + x[1] + x[2],
        ],
        [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
        domain_shape=(3,),
    )
    assert val.shape == (3,)


if __name__ == "__main__":
    test_multidim()
    # scheme_ = quadpy.c1.Fejer2(20)
    # # scheme_ = quadpy.c1.Midpoint()
    # test_scheme(scheme_)
    # test_show(scheme_)
    # # import matplotlib.pyplot as plt
    # # plt.savefig('demo.png', transparent=True)
