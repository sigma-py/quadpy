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
    a = 1.0
    b = 2.0
    schemes = [
        quadrature.line.Midpoint(),
        ]
    for scheme in schemes:
        yield check_triangle_scheme, scheme, a, b


def check_triangle_scheme(scheme, a, b):
    f = _create_test_polynomial(degree=scheme.degree)
    exact_val = _integrate_exact(f, a, b)
    val = quadrature.line.integrate(f, a, b, scheme)
    numpy.testing.assert_allclose(val, exact_val)
    return
