import math
from distutils.version import LooseVersion

import numpy
import orthopy
import pytest
import scipy
import sympy
from mpmath import mp

import quadpy


def test_gauss_sympy():
    n = 3
    a = 0
    b = 0
    _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.jacobi(
        n, a, b, "monic", symbolic=True
    )
    points, weights = quadpy.tools.scheme_from_rc(alpha, beta, "sympy")

    assert points == [-sympy.sqrt(sympy.S(3) / 5), 0, +sympy.sqrt(sympy.S(3) / 5)]
    assert weights == [sympy.S(5) / 9, sympy.S(8) / 9, sympy.S(5) / 9]


def test_gauss_mpmath():
    n = 5
    a = 0
    b = 0
    _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.jacobi(
        n, a, b, "monic", symbolic=True
    )
    mp.dps = 50
    points, weights = quadpy.tools.scheme_from_rc(alpha, beta, "mpmath")

    tol = 1.0e-50
    s = mp.sqrt(5 + 2 * mp.sqrt(mp.mpf(10) / mp.mpf(7))) / 3
    t = mp.sqrt(5 - 2 * mp.sqrt(mp.mpf(10) / mp.mpf(7))) / 3
    assert (abs(points - [-s, -t, 0.0, +t, +s]) < tol).all()

    u = mp.mpf(128) / mp.mpf(225)
    v = (322 + 13 * mp.sqrt(70)) / 900
    w = (322 - 13 * mp.sqrt(70)) / 900
    assert (abs(weights - [w, v, u, v, w]) < tol).all()


def test_gauss_numpy():
    n = 5
    tol = 1.0e-14
    _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.legendre(
        n, "monic", symbolic=False
    )
    flt = numpy.vectorize(float)
    alpha = flt(alpha)
    beta = flt(beta)
    points, weights = quadpy.tools.scheme_from_rc(alpha, beta, "numpy")

    s = math.sqrt(5.0 + 2 * math.sqrt(10.0 / 7.0)) / 3.0
    t = math.sqrt(5.0 - 2 * math.sqrt(10.0 / 7.0)) / 3.0
    assert (abs(points - [-s, -t, 0.0, +t, +s]) < tol).all()

    u = 128.0 / 225.0
    v = (322.0 + 13 * math.sqrt(70)) / 900.0
    w = (322.0 - 13 * math.sqrt(70)) / 900.0
    assert (abs(weights - [w, v, u, v, w]) < tol).all()


@pytest.mark.skipif(
    LooseVersion(scipy.__version__) < LooseVersion("1.0.0"), reason="Requires SciPy 1.0"
)
def test_jacobi_reconstruction(tol=1.0e-14):
    _, _, alpha1, beta1 = orthopy.line_segment.recurrence_coefficients.jacobi(
        4, 2, 1, "monic", symbolic=False
    )
    points, weights = quadpy.tools.scheme_from_rc(alpha1, beta1, "numpy")

    alpha2, beta2 = quadpy.tools.coefficients_from_gauss(points, weights)

    assert numpy.all(abs(alpha1 - alpha2) < tol)
    assert numpy.all(abs(beta1 - beta2) < tol)


@pytest.mark.skipif(
    LooseVersion(scipy.__version__) < LooseVersion("1.0.0"), reason="Requires SciPy 1.0"
)
def test_gautschi_how_to_and_how_not_to():
    """Test Gautschi's famous example from

    W. Gautschi,
    How and how not to check Gaussian quadrature formulae,
    BIT Numerical Mathematics,
    June 1983, Volume 23, Issue 2, pp 209–216,
    <https://doi.org/10.1007/BF02218441>.
    """
    points = numpy.array(
        [
            1.457697817613696e-02,
            8.102669876765460e-02,
            2.081434595902250e-01,
            3.944841255669402e-01,
            6.315647839882239e-01,
            9.076033998613676e-01,
            1.210676808760832,
            1.530983977242980,
            1.861844587312434,
            2.199712165681546,
            2.543839804028289,
            2.896173043105410,
            3.262066731177372,
            3.653371887506584,
            4.102376773975577,
        ]
    )
    weights = numpy.array(
        [
            3.805398607861561e-2,
            9.622028412880550e-2,
            1.572176160500219e-1,
            2.091895332583340e-1,
            2.377990401332924e-1,
            2.271382574940649e-1,
            1.732845807252921e-1,
            9.869554247686019e-2,
            3.893631493517167e-2,
            9.812496327697071e-3,
            1.439191418328875e-3,
            1.088910025516801e-4,
            3.546866719463253e-6,
            3.590718819809800e-8,
            5.112611678291437e-11,
        ]
    )

    # weight function exp(-t**3/3)
    n = len(points)
    moments = numpy.array(
        [3.0 ** ((k - 2) / 3.0) * math.gamma((k + 1) / 3.0) for k in range(2 * n)]
    )

    alpha, beta = quadpy.tools.coefficients_from_gauss(points, weights)
    # alpha, beta = quadpy.tools.chebyshev(moments)

    errors_alpha, errors_beta = orthopy.tools.gautschi_test_3(moments, alpha, beta)

    assert numpy.max(errors_alpha) > 1.0e-2
    assert numpy.max(errors_beta) > 1.0e-2


# def test_expt3():
#     '''Full example from Gautschi's "How to and how not to" article.
#     '''
#     # moments = quadpy.tools.integrate(
#     #         lambda x: sympy.exp(-x**3/3),
#     #         0, sympy.oo,
#     #         31
#     #         )
#     # print(moments)
#     # alpha, beta = quadpy.tools.chebyshev(moments)
#
#     alpha, beta = quadpy.tools.stieltjes(
#             lambda x: sympy.exp(-x**3/3),
#             0, sympy.oo,
#             5
#             )
#     print(alpha)
#     print(beta)


@pytest.mark.parametrize("k", [0, 2, 4])
def test_xk(k):
    n = 10

    moments = quadpy.tools.integrate(
        lambda x: [x ** (i + k) for i in range(2 * n)], -1, +1
    )

    alpha, beta = quadpy.tools.chebyshev(moments)

    assert (alpha == 0).all()
    assert beta[0] == moments[0]
    assert beta[1] == sympy.S(k + 1) / (k + 3)
    assert beta[2] == sympy.S(4) / ((k + 5) * (k + 3))
    quadpy.tools.scheme_from_rc(
        numpy.array([sympy.N(a) for a in alpha], dtype=float),
        numpy.array([sympy.N(b) for b in beta], dtype=float),
        mode="numpy",
    )

    def leg_polys(x):
        return orthopy.line_segment.tree_legendre(x, 19, "monic", symbolic=True)

    moments = quadpy.tools.integrate(
        lambda x: [x ** k * leg_poly for leg_poly in leg_polys(x)], -1, +1
    )

    _, _, a, b = orthopy.line_segment.recurrence_coefficients.legendre(
        2 * n, "monic", symbolic=True
    )

    alpha, beta = quadpy.tools.chebyshev_modified(moments, a, b)

    assert (alpha == 0).all()
    assert beta[0] == moments[0]
    assert beta[1] == sympy.S(k + 1) / (k + 3)
    assert beta[2] == sympy.S(4) / ((k + 5) * (k + 3))
    points, weights = quadpy.tools.scheme_from_rc(
        numpy.array([sympy.N(a) for a in alpha], dtype=float),
        numpy.array([sympy.N(b) for b in beta], dtype=float),
        mode="numpy",
    )


if __name__ == "__main__":
    # test_gauss('mpmath')
    # test_logo()
    test_xk(2)
