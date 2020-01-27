import math

import numpy
import pytest
from mpmath import mp

import quadpy
from helpers import check_degree_1d


@pytest.mark.parametrize(
    "scheme",
    [quadpy.line_segment.midpoint()]
    + [quadpy.line_segment.trapezoidal()]
    + [quadpy.line_segment.clenshaw_curtis(k) for k in range(2, 10)]
    + [quadpy.line_segment.gauss_legendre(k) for k in range(1, 6)]
    + [quadpy.line_segment.gauss_lobatto(k) for k in range(2, 7)]
    + [quadpy.line_segment.gauss_kronrod(k) for k in range(2, 7)]
    + [quadpy.line_segment.gauss_patterson(k) for k in range(9)]
    + [quadpy.line_segment.gauss_radau(k) for k in range(2, 10)]
    + [quadpy.line_segment.fejer_1(k) for k in range(1, 10)]
    + [quadpy.line_segment.fejer_2(k) for k in range(1, 10)]
    + [quadpy.line_segment.newton_cotes_closed(k) for k in range(1, 5)]
    + [quadpy.line_segment.newton_cotes_open(k) for k in range(1, 5)],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    # https://github.com/nschloe/quadpy/issues/227
    assert scheme.weights.ndim == 1

    degree = 0
    while True:
        # Set bounds such that the values are between 0.5 and 1.5.
        exact_val = 1.0 / (degree + 1)
        interval = numpy.array(
            [
                [0.5 ** (1.0 / (degree + 1)), 0.0, 0.0],
                [1.5 ** (1.0 / (degree + 1)), 0.0, 0.0],
            ]
        )
        interval = numpy.array([[0.3], [0.5]])
        val = scheme.integrate(lambda x: x[0] ** degree, interval)
        # same test with line embedded in R^2
        interval = numpy.array(
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
    "scheme", [quadpy.line_segment.chebyshev_gauss_1(k) for k in range(1, 10)]
)
def test_cheb1_scheme(scheme):
    def integrate_exact(k):
        # \int_-1^1 x^k / sqrt(1 - x^2)
        if k == 0:
            return numpy.pi
        if k % 2 == 1:
            return 0.0
        return (
            numpy.sqrt(numpy.pi)
            * ((-1) ** k + 1)
            * math.gamma(0.5 * (k + 1))
            / math.gamma(0.5 * k)
            / k
        )

    degree = check_degree_1d(
        lambda poly: scheme.integrate(poly, numpy.array([[-1.0], [1.0]])),
        integrate_exact,
        scheme.degree + 1,
    )
    assert degree >= scheme.degree


@pytest.mark.parametrize(
    "scheme", [quadpy.line_segment.chebyshev_gauss_2(k) for k in range(1, 10)]
)
def test_cheb2_scheme(scheme):
    def integrate_exact(k):
        # \int_-1^1 x^k * sqrt(1 - x^2)
        if k == 0:
            return 0.5 * numpy.pi
        if k % 2 == 1:
            return 0.0
        return (
            numpy.sqrt(numpy.pi)
            * ((-1) ** k + 1)
            * math.gamma(0.5 * (k + 1))
            / math.gamma(0.5 * k + 2)
            / 4
        )

    degree = check_degree_1d(
        lambda poly: scheme.integrate(poly, numpy.array([[-1.0], [1.0]])),
        integrate_exact,
        scheme.degree + 1,
    )
    assert degree >= scheme.degree


@pytest.mark.parametrize("scheme", [quadpy.line_segment.newton_cotes_closed(5)])
def test_show(scheme):
    scheme.show()


def test_integrate_split():
    val = quadpy.line_segment.trapezoidal().integrate_split(
        lambda r: 0.5108
        / r ** 2
        / numpy.sqrt(2 * 1.158 + 2 / r - 0.5108 ** 2 / (2 * r ** 2)),
        0.15,
        0.702,
        100,
    )
    reference = 0.961715
    assert abs(val - reference) < 1.0e-3 * reference


def test_legendre_mpmath():
    mp.dps = 50
    scheme = quadpy.line_segment.gauss_legendre(4, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.sqrt(mp.mpf(3) / 7 - mp.mpf(2) / 7 * mp.sqrt(mp.mpf(6) / 5))
    x2 = mp.sqrt(mp.mpf(3) / 7 + mp.mpf(2) / 7 * mp.sqrt(mp.mpf(6) / 5))
    assert (abs(scheme.points - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = (18 + mp.sqrt(30)) / 36
    w2 = (18 - mp.sqrt(30)) / 36
    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()


def test_chebyshev1_sympy():
    scheme = quadpy.line_segment.chebyshev_gauss_1(4, mode="sympy")
    scheme_numpy = quadpy.line_segment.chebyshev_gauss_1(4, mode="numpy")

    flt = numpy.vectorize(float)
    tol = 1.0e-15

    assert (abs(flt(scheme.points) - scheme_numpy.points) < tol).all()
    assert (abs(flt(scheme.weights) - scheme_numpy.weights) < tol).all()


def test_chebyshev2_sympy():
    scheme = quadpy.line_segment.chebyshev_gauss_2(4, mode="sympy")
    scheme_numpy = quadpy.line_segment.chebyshev_gauss_2(4, mode="numpy")

    flt = numpy.vectorize(float)
    tol = 1.0e-15

    assert (abs(flt(scheme.points) - scheme_numpy.points) < tol).all()
    assert (abs(flt(scheme.weights) - scheme_numpy.weights) < tol).all()


def test_chebyshev1_mpmath():
    mp.dps = 50
    scheme = quadpy.line_segment.chebyshev_gauss_1(4, mode="mpmath")
    tol = 1.0e-50

    x1 = mp.cos(3 * mp.pi / 8)
    x2 = mp.cos(1 * mp.pi / 8)
    assert (abs(scheme.points - [+x2, +x1, -x1, -x2]) < tol).all()

    w = mp.pi / 4
    tol = 1.0e-49
    assert (abs(scheme.weights - [w, w, w, w]) < tol).all()


def test_chebyshev2_mpmath():
    mp.dps = 51
    scheme = quadpy.line_segment.chebyshev_gauss_2(4, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.cos(2 * mp.pi / 5)
    x2 = mp.cos(1 * mp.pi / 5)
    assert (abs(scheme.points - [+x2, +x1, -x1, -x2]) < tol).all()

    w1 = mp.pi / 5 * mp.sin(2 * mp.pi / 5) ** 2
    w2 = mp.pi / 5 * mp.sin(1 * mp.pi / 5) ** 2
    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()


def test_jacobi_mpmath():
    mp.dps = 51
    scheme = quadpy.line_segment.gauss_jacobi(4, 1, 1, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.sqrt((7 - 2 * mp.sqrt(7)) / 21)
    x2 = mp.sqrt((7 + 2 * mp.sqrt(7)) / 21)
    assert (abs(scheme.points - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = (5 + mp.sqrt(7)) / 15
    w2 = (5 - mp.sqrt(7)) / 15
    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()


def test_multidim():
    scheme = quadpy.line_segment.gauss_legendre(5)

    # simple scalar integration
    val = scheme.integrate(numpy.sin, [0.0, 1.0])
    assert val.shape == ()

    # scalar integration on 3 subdomains
    val = scheme.integrate(numpy.sin, [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]])
    assert val.shape == (3,)

    # scalar integration in 3D
    val = scheme.integrate(
        lambda x: x[0] + numpy.sin(x[1]) + numpy.cos(x[2]),
        [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
    )
    assert val.shape == ()

    # vector-valued integration on 3 subdomains
    val = scheme.integrate(
        lambda x: [numpy.sin(x), numpy.cos(x)],
        [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]]
    )
    assert val.shape == (2, 3)

    # vector-valued integration in 3D
    val = scheme.integrate(
        lambda x: [x[0] + numpy.sin(x[1]), numpy.cos(x[0]) * x[2]],
        [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]]
    )
    assert val.shape == (2,)


if __name__ == "__main__":
    test_multidim()
    # scheme_ = quadpy.line_segment.Fejer2(20)
    # # scheme_ = quadpy.line_segment.Midpoint()
    # test_scheme(scheme_)
    # test_show(scheme_)
    # # import matplotlib.pyplot as plt
    # # plt.savefig('demo.png', transparent=True)
