# -*- coding: utf-8 -*-
#
import numpy
import numpy.testing
import quadrature
import sympy

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def test_generator():
    schemes = [
        quadrature.line.Midpoint(),
        quadrature.line.Trapezoidal(),
        quadrature.line.GaussLegendre(1),
        quadrature.line.GaussLegendre(2),
        quadrature.line.GaussLegendre(3),
        quadrature.line.GaussLegendre(4),
        quadrature.line.GaussLegendre(5),
        quadrature.line.GaussPatterson(0),
        quadrature.line.GaussPatterson(1),
        quadrature.line.GaussPatterson(2),
        quadrature.line.GaussPatterson(3),
        quadrature.line.GaussPatterson(4),
        quadrature.line.GaussPatterson(5),
        quadrature.line.GaussPatterson(6),
        quadrature.line.ClenshawCurtis(1),
        quadrature.line.ClenshawCurtis(2),
        quadrature.line.ClenshawCurtis(3),
        quadrature.line.ClenshawCurtis(4),
        quadrature.line.ClenshawCurtis(5),
        quadrature.line.ClenshawCurtis(9),
        quadrature.line.ClenshawCurtis(17),
        quadrature.line.ClenshawCurtis(33),
        quadrature.line.ClenshawCurtis(65),
        quadrature.line.NewtonCotesClosed(1),
        quadrature.line.NewtonCotesClosed(2),
        quadrature.line.NewtonCotesClosed(3),
        quadrature.line.NewtonCotesClosed(4),
        quadrature.line.NewtonCotesOpen(2),
        quadrature.line.NewtonCotesOpen(3),
        quadrature.line.NewtonCotesOpen(4),
        quadrature.line.NewtonCotesOpen(5),
        ]
    for scheme in schemes:
        yield check_scheme, scheme


def check_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    degree = 0
    while True:
        def poly(x): return x**degree
        # Set bounds such that the values are between 0.5 and 1.5.
        a = 0.5**(1.0/(degree+1))
        b = 1.5**(1.0/(degree+1))
        exact_val = 1.0/(degree+1)
        val = quadrature.line.integrate(poly, a, b, scheme)
        if abs(exact_val - val) / abs(exact_val) > 1.0e-12:
            break
        if degree >= scheme.degree:
            break
        degree += 1
    assert degree >= scheme.degree
    return


def test_show():
    quadrature.line.show(
        0.0, 1.0,
        # quadrature.line.NewtonCotesOpen(4)
        # quadrature.line.GaussLegendre(31)
        # quadrature.line.GaussPatterson(4)
        quadrature.line.ClenshawCurtis(33)
        )
    return


if __name__ == '__main__':
    # test_show()
    # plt.show()
    scheme = quadrature.line.GaussLegendre(7)
    check_scheme(scheme)
