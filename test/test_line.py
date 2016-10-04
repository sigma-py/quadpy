# -*- coding: utf-8 -*-
#
import numpy
import numpy.testing
import quadrature
import sympy


def _integrate_exact(f, a, b):
    x = sympy.Symbol('x')
    exact = sympy.integrate(f(x), (x, a, b))
    return float(exact)


def _create_test_polynomial(degree):
    '''Return polynomials of the form

      p0(x) = 1,
      p1(x) = 1 + 1/2 * x
      p2(x) = 1 + 1/2 * x + 1/3 * x**2

    etc.
    '''
    def f(x):
        return numpy.sum([1.0/(k+1) * x**k for k in range(degree+1)])
    return f


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
        ]
    for scheme in schemes:
        yield check_triangle_scheme, scheme, a, b


def check_triangle_scheme(scheme, a, b):
    f = _create_test_polynomial(degree=scheme.degree)
    exact_val = _integrate_exact(f, a, b)
    val = quadrature.line.integrate(f, a, b, scheme)
    numpy.testing.assert_allclose(val, exact_val)
    return
