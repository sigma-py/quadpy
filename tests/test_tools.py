import math

import numpy as np
import orthopy
import pytest
import scipy
import sympy
from mpmath import mp
from packaging import version

import quadpy


def test_gauss_sympy():
    n = 3
    rc = orthopy.c1.jacobi.RecurrenceCoefficients(
        "monic", alpha=0, beta=0, symbolic=True
    )
    _, a, b = np.array([rc[k] for k in range(n)]).T
    points, weights = quadpy.tools.scheme_from_rc(a, b, rc.int_1, "sympy")

    assert points == [-sympy.sqrt(sympy.S(3) / 5), 0, +sympy.sqrt(sympy.S(3) / 5)]
    assert weights == [sympy.S(5) / 9, sympy.S(8) / 9, sympy.S(5) / 9]


def test_gauss_mpmath():
    n = 5

    rc = orthopy.c1.jacobi.RecurrenceCoefficients(
        "monic", alpha=0, beta=0, symbolic=True
    )
    _, alpha, beta = np.array([rc[k] for k in range(n)]).T
    mp.dps = 50
    points, weights = quadpy.tools.scheme_from_rc(alpha, beta, rc.int_1, "mpmath")

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
    rc = orthopy.c1.legendre.RecurrenceCoefficients("monic", symbolic=False)
    _, alpha, beta = np.array([rc[k] for k in range(n)]).T

    flt = np.vectorize(float)
    alpha = flt(alpha)
    beta = flt(beta)
    points, weights = quadpy.tools.scheme_from_rc(alpha, beta, rc.int_1, "numpy")

    s = math.sqrt(5.0 + 2 * math.sqrt(10.0 / 7.0)) / 3.0
    t = math.sqrt(5.0 - 2 * math.sqrt(10.0 / 7.0)) / 3.0
    assert (abs(points - [-s, -t, 0.0, +t, +s]) < tol).all()

    u = 128.0 / 225.0
    v = (322.0 + 13 * math.sqrt(70)) / 900.0
    w = (322.0 - 13 * math.sqrt(70)) / 900.0
    assert (abs(weights - [w, v, u, v, w]) < tol).all()


@pytest.mark.skipif(
    version.parse(scipy.__version__) < version.parse("1.0.0"),
    reason="requires scipy 1.0 or higher",
)
def test_jacobi_reconstruction(tol=1.0e-14):
    n = 4
    rc = orthopy.c1.jacobi.RecurrenceCoefficients(
        "monic", alpha=2, beta=1, symbolic=False
    )
    _, alpha1, beta1 = np.array([rc[k] for k in range(n)]).T

    points, weights = quadpy.tools.scheme_from_rc(alpha1, beta1, rc.int_1, "numpy")

    alpha2, beta2, int_1 = quadpy.tools.coefficients_from_gauss(points, weights)

    assert np.all(abs(alpha1 - alpha2) < tol)
    assert np.all(abs(beta1[1:] - beta2[1:]) < tol)
    assert abs(rc.int_1 == int_1) < tol


@pytest.mark.skip(reason="wait for new orthopy")
def test_gautschi_how_to_and_how_not_to():
    """Test Gautschi's famous example from

    W. Gautschi,
    How and how not to check Gaussian quadrature formulae,
    BIT Numerical Mathematics,
    June 1983, Volume 23, Issue 2, pp 209–216,
    <https://doi.org/10.1007/BF02218441>.
    """
    points = np.array(
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
    weights = np.array(
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
    moments = np.array(
        [3.0 ** ((k - 2) / 3.0) * math.gamma((k + 1) / 3.0) for k in range(2 * n)]
    )

    alpha, beta, _ = quadpy.tools.coefficients_from_gauss(points, weights)
    # alpha, beta = quadpy.tools.chebyshev(moments)

    errors_alpha, errors_beta = orthopy.tools.gautschi_test_3(moments, alpha, beta)

    assert np.max(errors_alpha) > 1.0e-2
    assert np.max(errors_beta) > 1.0e-2


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

    x = sympy.Symbol("x")
    moments = [sympy.integrate(x ** (i + k), (x, -1, 1)) for i in range(2 * n)]

    alpha, beta, int_1 = orthopy.tools.chebyshev(moments)

    assert (alpha == 0).all()
    assert math.isnan(beta[0])
    assert int_1 == moments[0]
    assert beta[1] == sympy.S(k + 1) / (k + 3)
    assert beta[2] == sympy.S(4) / ((k + 5) * (k + 3))
    quadpy.tools.scheme_from_rc(
        np.array([sympy.N(a) for a in alpha], dtype=float),
        np.array([sympy.N(b) for b in beta], dtype=float),
        int_1,
        mode="numpy",
    )

    evaluator = orthopy.c1.legendre.Eval(x, "monic", symbolic=True)
    moments = [
        sympy.integrate(x ** k * next(evaluator), (x, -1, 1)) for _ in range(2 * n)
    ]
    rc = orthopy.c1.legendre.RecurrenceCoefficients("monic", symbolic=True)

    alpha, beta, int_1 = orthopy.tools.chebyshev_modified(moments, rc)

    assert (alpha == 0).all()
    assert math.isnan(beta[0])
    assert int_1 == moments[0]
    assert beta[1] == sympy.S(k + 1) / (k + 3)
    assert beta[2] == sympy.S(4) / ((k + 5) * (k + 3))
    points, weights = quadpy.tools.scheme_from_rc(
        np.array([sympy.N(a) for a in alpha], dtype=float),
        np.array([sympy.N(b) for b in beta], dtype=float),
        int_1,
        mode="numpy",
    )


if __name__ == "__main__":
    # test_gauss('mpmath')
    # test_logo()
    test_xk(2)
