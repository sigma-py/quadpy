import math
from collections import namedtuple

import numpy as np
import orthopy

from ..helpers import article
from ..tools import scheme_from_rc
from ._gauss_legendre import gauss_legendre
from ._helpers import C1Scheme, _find_shapes

source = article(
    authors=["Dirk P. Laurie"],
    title="Calculation of Gauss-Kronrod quadrature rules",
    journal="Math. Comp.",
    volume="66",
    year="1997",
    pages="1133-1145",
    url="https://doi.org/10.1090/S0025-5718-97-00861-2",
)


def gauss_kronrod(n, a=0, b=0):
    """
    Gauss-Kronrod quadrature; see
    <https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula>.

    Besides points and weights, this class provides the weights of the corresponding
    Gauss-Legendre scheme in self.gauss_weights.

    Code adapted from
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/r_kronrod.m>,
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/kronrod.m>.
    """
    # The general scheme is:
    # Get the Jacobi recurrence coefficients, get the Kronrod vectors alpha and beta,
    # and hand those off to scheme_from_rc. There, the eigenproblem for a tridiagonal
    # matrix with alpha and beta is solved to retrieve the points and weights.
    # TODO replace math.ceil by -(-k//n)
    length = int(math.ceil(3 * n / 2.0)) + 1
    degree = 2 * length + 1

    rc = orthopy.c1.jacobi.RecurrenceCoefficients("monic", a, b, symbolic=False)
    _, alpha, beta = np.array([rc[k] for k in range(length)]).T
    beta[0] = rc.int_1

    a, b = _r_kronrod(n, alpha, beta)
    x, w = scheme_from_rc(a, b, rc.int_1, mode="numpy")
    # sort by x
    i = np.argsort(x)
    points = x[i]
    weights = w[i]
    return C1Scheme(f"Gauss-Kronrod ({n})", degree, weights, points, source)


def _r_kronrod(n, a0, b0):
    assert len(a0) == int(math.ceil(3 * n / 2.0)) + 1
    assert len(b0) == int(math.ceil(3 * n / 2.0)) + 1

    a = np.zeros(2 * n + 1)
    b = np.zeros(2 * n + 1)

    k = int(math.floor(3 * n / 2.0)) + 1
    a[:k] = a0[:k]
    k = int(math.ceil(3 * n / 2.0)) + 1
    b[:k] = b0[:k]
    s = np.zeros(int(math.floor(n / 2.0)) + 2)
    t = np.zeros(int(math.floor(n / 2.0)) + 2)
    t[1] = b[n + 1]
    for m in range(n - 1):
        k0 = int(math.floor((m + 1) / 2.0))
        k = np.arange(k0, -1, -1)
        L = m - k
        s[k + 1] = np.cumsum(
            (a[k + n + 1] - a[L]) * t[k + 1] + b[k + n + 1] * s[k] - b[L] * s[k + 1]
        )
        s, t = t, s

    j = int(math.floor(n / 2.0)) + 1
    s[1 : j + 1] = s[:j]
    for m in range(n - 1, 2 * n - 2):
        k0 = m + 1 - n
        k1 = int(math.floor((m - 1) / 2.0))
        k = np.arange(k0, k1 + 1)
        L = m - k
        j = n - 1 - L
        s[j + 1] = np.cumsum(
            -(a[k + n + 1] - a[L]) * t[j + 1]
            - b[k + n + 1] * s[j + 1]
            + b[L] * s[j + 2]
        )
        j = j[-1]
        k = int(math.floor((m + 1) / 2.0))
        if m % 2 == 0:
            a[k + n + 1] = a[k] + (s[j + 1] - b[k + n + 1] * s[j + 2]) / t[j + 2]
        else:
            b[k + n + 1] = s[j + 1] / s[j + 2]
        s, t = t, s

    a[2 * n] = a[n - 1] - b[2 * n] * s[1] / t[1]
    return a, b


def _gauss_kronrod_integrate(
    k,
    f,
    intervals,
    dot=np.dot,
    domain_shape=None,
    range_shape=None,
):
    # Compute the integral estimations according to Gauss and Gauss-Kronrod, sharing the
    # function evaluations
    gk = gauss_kronrod(k)
    gl = gauss_legendre(k)
    # scale points
    x0 = 0.5 * (1.0 - gk.points)
    x1 = 0.5 * (1.0 + gk.points)

    sp = np.multiply.outer(intervals[0], x0) + np.multiply.outer(intervals[1], x1)
    if domain_shape is not None and range_shape is not None:
        # Make it easy on the f by flattening out interval_set
        interval_set_shape = intervals.shape[1 + len(domain_shape) :]
        sp = sp.reshape(*domain_shape, -1)
        fx_gk = np.asarray(f(sp))
        fx_gk = fx_gk.reshape(*range_shape, *interval_set_shape, *x0.shape)
    else:
        fx_gk = np.asarray(f(sp))
        # try and guess shapes of domain, range, intervals
        domain_shape, range_shape, interval_set_shape = _find_shapes(
            fx_gk, intervals, gk.points, domain_shape, range_shape
        )

    fx_gl = fx_gk[..., 1::2]

    diff = intervals[1] - intervals[0]
    interval_lengths = np.sqrt(np.sum(diff ** 2, axis=tuple(range(len(domain_shape)))))
    assert interval_lengths.shape == interval_set_shape

    # integrate
    # average value of gk over the interval
    average_gk = 0.5 * dot(fx_gk, gk.weights)
    val_gauss_kronrod = average_gk * interval_lengths
    val_gauss_legendr = 0.5 * interval_lengths * dot(fx_gl, gl.weights)

    assert val_gauss_kronrod.shape == range_shape + interval_set_shape
    assert val_gauss_legendr.shape == range_shape + interval_set_shape

    # Get an error estimate. According to
    #
    #   Pedro Gonnet,
    #   A Review of Error Estimation in Adaptive Quadrature,
    #   ACM Computing Surveys (CSUR) Surveys,
    #   Volume 44, Issue 4, August 2012
    #   <https://doi.org/10.1145/2333112.2333117>,
    #   <https://arxiv.org/pdf/1003.4629.pdf>
    #
    # the classical QUADPACK still compares favorably with other approaches.
    fx_avg_abs = np.abs(fx_gk - average_gk[..., None])
    I_tilde = 0.5 * interval_lengths * dot(fx_avg_abs, gk.weights)

    # The exponent 1.5 is chosen such that (200*x)**1.5 is approximately x at 1.0e-6,
    # the machine precision on IEEE 754 32-bit floating point arithmentic. This could be
    # adapted to
    #
    #   eps = np.finfo(float).eps
    #   exponent = np.log(eps) / np.log(200*eps)
    #
    # The following expression is
    #
    # error_estimate = I_tilde * np.minimum(
    #     np.ones(I_tilde.shape),
    #     (200 * abs(val_gauss_kronrod - val_gauss_legendr) / I_tilde) ** 1.5,
    # )
    #
    # with handling NaNs (if I_tilde is 0).
    error_estimate = np.empty(I_tilde.shape)
    idx = np.abs(I_tilde) > 1.0e-15
    vals = (
        200 * np.abs(val_gauss_kronrod[idx] - val_gauss_legendr[idx]) / I_tilde[idx]
    ) ** 1.5
    error_estimate[idx] = np.minimum(np.ones(np.sum(idx)), vals)
    error_estimate[~idx] = 1.0
    error_estimate *= I_tilde

    assert error_estimate.shape == range_shape + interval_set_shape

    gk_return = namedtuple(
        "gauss_kronrod_return_values",
        [
            "val_gauss_kronrod",
            "val_gauss_legendre",
            "interval_lengths",
            "error_estimate",
            "domain_shape",
            "range_shape",
        ],
    )

    return gk_return(
        val_gauss_kronrod,
        val_gauss_legendr,
        interval_lengths,
        error_estimate,
        domain_shape,
        range_shape,
    )
