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


def _integrate_exact(f, a, b):
    x = sympy.Symbol('x')
    exact = sympy.integrate(f(x), (x, a, b))
    return float(exact)


def test_generator():
    a = -0.3
    b = 0.5
    schemes = [
        quadrature.line.Midpoint(),
        quadrature.line.Trapezoidal(),
        quadrature.line.GaussLegendre(1),
        quadrature.line.GaussLegendre(2),
        quadrature.line.GaussLegendre(3),
        quadrature.line.GaussLegendre(4),
        quadrature.line.GaussLegendre(7),
        quadrature.line.GaussLegendre(8),
        quadrature.line.GaussLegendre(15),
        quadrature.line.GaussLegendre(16),
        quadrature.line.GaussLegendre(31),
        quadrature.line.GaussLegendre(32),
        quadrature.line.GaussLegendre(63),
        quadrature.line.GaussLegendre(64),
        quadrature.line.GaussLegendre(65),
        quadrature.line.GaussLegendre(127),
        quadrature.line.GaussPatterson(1),
        quadrature.line.GaussPatterson(3),
        quadrature.line.GaussPatterson(7),
        quadrature.line.GaussPatterson(15),
        quadrature.line.GaussPatterson(31),
        quadrature.line.GaussPatterson(63),
        quadrature.line.GaussPatterson(127),
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
        yield check_scheme, scheme, a, b


def check_scheme(scheme, a, b):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    degree = 0
    while True:
        def poly(x): return x**degree
        exact_val = _integrate_exact(poly, a, b)
        val = quadrature.line.integrate(poly, a, b, scheme)
        if abs(exact_val - val) > 1.0e-10:
            break
        degree += 1
    numpy.testing.assert_equal(degree-1, scheme.degree)
    return


def test_show():
    quadrature.line.show(
        0.0, 1.0,
        # quadrature.line.NewtonCotesOpen(4)
        # quadrature.line.GaussLegendre(31)
        # quadrature.line.GaussPatterson(31)
        quadrature.line.ClenshawCurtis(33)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
