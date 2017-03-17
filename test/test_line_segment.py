# -*- coding: utf-8 -*-
#
from helpers import check_degree_1d

import math
from matplotlib import pyplot as plt
import numpy
import pytest
import quadpy


@pytest.mark.parametrize(
    'scheme',
    [quadpy.line_segment.Midpoint()]
    + [quadpy.line_segment.Trapezoidal()]
    + [quadpy.line_segment.ClenshawCurtis(k) for k in range(2, 10)]
    + [quadpy.line_segment.GaussLegendre(k) for k in range(1, 6)]
    + [quadpy.line_segment.GaussLobatto(k) for k in range(2, 7)]
    + [quadpy.line_segment.GaussKronrod(k) for k in range(2, 7)]
    + [quadpy.line_segment.GaussPatterson(k) for k in range(7)]
    + [quadpy.line_segment.GaussRadau(k) for k in range(2, 10)]
    + [quadpy.line_segment.Fejer1(k) for k in range(1, 10)]
    + [quadpy.line_segment.Fejer2(k) for k in range(1, 10)]
    + [quadpy.line_segment.NewtonCotesClosed(k) for k in range(1, 5)]
    + [quadpy.line_segment.NewtonCotesOpen(k) for k in range(1, 5)]
    )
def test_scheme(scheme):
    degree = 0
    while True:
        # Set bounds such that the values are between 0.5 and 1.5.
        interval = [
            0.5**(1.0/(degree+1)),
            1.5**(1.0/(degree+1)),
            ]
        exact_val = 1.0/(degree+1)
        val = quadpy.line_segment.integrate(
                lambda x: x**degree,
                interval, scheme
                )
        if abs(exact_val - val) / abs(exact_val) > 1.0e-12:
            break
        if degree >= scheme.degree:
            break
        degree += 1
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.line_segment.ChebyshevGauss1(k) for k in range(1, 10)]
    )
def test_cheb1_scheme(scheme):
    def integrate_exact(k):
        # \int_-1^1 x^k / sqrt(1 - x^2)
        if k == 0:
            return numpy.pi
        if k % 2 == 1:
            return 0.0
        return numpy.sqrt(numpy.pi) * ((-1)**k + 1) \
            * math.gamma(0.5*(k+1)) / math.gamma(0.5*k) \
            / k

    degree = check_degree_1d(
            lambda poly: quadpy.line_segment.integrate(
                    poly, [-1.0, 1.0], scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.line_segment.ChebyshevGauss2(k) for k in range(1, 10)]
    )
def test_cheb2_scheme(scheme):
    def integrate_exact(k):
        # \int_-1^1 x^k * sqrt(1 - x^2)
        if k == 0:
            return 0.5 * numpy.pi
        if k % 2 == 1:
            return 0.0
        return numpy.sqrt(numpy.pi) * ((-1)**k + 1) \
            * math.gamma(0.5*(k+1)) / math.gamma(0.5*k+2) \
            / 4

    degree = check_degree_1d(
            lambda poly: quadpy.line_segment.integrate(
                    poly, [-1.0, 1.0], scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.line_segment.GaussLaguerre(k) for k in range(1, 10)]
    )
def test_laguerre_scheme(scheme):
    def integrate_exact(k):
        # \int_0^\infty x^k * exp(-x)
        return math.gamma(k + 1)

    degree = check_degree_1d(
            lambda poly: quadpy.line_segment.integrate(
                    poly, [-1.0, 1.0], scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.line_segment.GaussHermite(k) for k in range(1, 8)]
    )
def test_hermite_scheme(scheme):
    def integrate_exact(k):
        # \int_-\infty^\infty x^k * exp(-x^2)
        return 0.5 * ((-1)**k + 1) * math.gamma(0.5*(k + 1))

    degree = check_degree_1d(
            lambda poly: quadpy.line_segment.integrate(
                    poly, [-1.0, 1.0], scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.line_segment.NewtonCotesClosed(5)]
    )
def test_show(scheme):
    quadpy.line_segment.show(scheme)
    return


if __name__ == '__main__':
    scheme = quadpy.line_segment.Fejer2(20)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
