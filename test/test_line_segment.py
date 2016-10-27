# -*- coding: utf-8 -*-
#
from helpers import check_degree_1d

import math
import numpy
import quadrature

import os
import matplotlib as mpl
import pytest
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.Midpoint()]
    + [quadrature.line_segment.Trapezoidal()]
    + [quadrature.line_segment.GaussLegendre(k) for k in range(1, 6)]
    + [quadrature.line_segment.GaussLobatto(k) for k in range(2, 7)]
    + [quadrature.line_segment.GaussPatterson(k) for k in range(7)]
    + [quadrature.line_segment.ClenshawCurtis(k) for k in [
        1, 2, 3, 4, 5, 9, 17, 33, 65
        ]]
    + [quadrature.line_segment.NewtonCotesClosed(k) for k in range(1, 5)]
    + [quadrature.line_segment.NewtonCotesOpen(k) for k in range(1, 5)]
    )
def test_scheme(scheme):
    degree = 0
    while True:
        # Set bounds such that the values are between 0.5 and 1.5.
        a = 0.5**(1.0/(degree+1))
        b = 1.5**(1.0/(degree+1))
        exact_val = 1.0/(degree+1)
        val = quadrature.line_segment.integrate(
                lambda x: x**degree,
                a, b, scheme
                )
        if abs(exact_val - val) / abs(exact_val) > 1.0e-12:
            break
        if degree >= scheme.degree:
            break
        degree += 1
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.ChebyshevGauss1(k) for k in range(1, 10)]
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
            lambda poly: quadrature.line_segment.integrate(
                    poly, -1.0, 1.0, scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.ChebyshevGauss2(k) for k in range(1, 10)]
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
            lambda poly: quadrature.line_segment.integrate(
                    poly, -1.0, 1.0, scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.GaussLaguerre(k) for k in range(1, 10)]
    )
def test_laguerre_scheme(scheme):
    def integrate_exact(k):
        # \int_0^\infty x^k * exp(-x)
        return math.gamma(k + 1)

    degree = check_degree_1d(
            lambda poly: quadrature.line_segment.integrate(
                    poly, -1.0, 1.0, scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.GaussHermite(k) for k in range(1, 10)]
    )
def test_hermite_scheme(scheme):
    def integrate_exact(k):
        # \int_-\infty^\infty x^k * exp(-x^2)
        return 0.5 * ((-1)**k + 1) * math.gamma(0.5*(k + 1))

    degree = check_degree_1d(
            lambda poly: quadrature.line_segment.integrate(
                    poly, -1.0, 1.0, scheme
                    ),
            integrate_exact,
            lambda degree: [[degree]],
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.NewtonCotesClosed(5)]
    )
def test_show(scheme):
    quadrature.line_segment.show(0.0, 1.0, scheme)
    return


if __name__ == '__main__':
    scheme = quadrature.line_segment.GaussLobatto(2)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
