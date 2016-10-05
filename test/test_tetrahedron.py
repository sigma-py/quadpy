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


def _integrate_exact(f, tetrahedron):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference tetrahedron [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the tetrahedron. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    def g(xi):
        pxi = tetrahedron[0] * (1 - xi[0] - xi[1] - xi[2]) \
            + tetrahedron[1] * xi[0] \
            + tetrahedron[2] * xi[1] \
            + tetrahedron[3] * xi[2]
        return f(pxi)

    x = sympy.DeferredVector('x')
    exact = 6 * quadrature.tetrahedron.volume(tetrahedron) \
        * sympy.integrate(
            sympy.integrate(
              sympy.integrate(g(x), (x[2], 0, 1-x[0]-x[1])),
              (x[1], 0, 1-x[0])
              ),
            (x[0], 0, 1)
          )
    return float(exact)


def _get_sum_tuples(a):
    '''Returns the list of 3-tuples with sum a.
    '''
    tuples = []
    for i in range(a+1):
        for j in range(a+1-i):
            k = a - i - j
            tuples.append((i, j, k))
    return tuples


def _create_test_polynomial(degree):
    '''The k-th order terms in polynomial have the form

        alpha_{k,i} * x^k0 * y^k1 * z^k2,

    with k0 + k1 + k2 = k. Take

      alpha_{k, i} = (i+1) / (k+2)

    such that

      p0(x) = 1/3,

      p1(x) = 1/3 \
            + 1/4 * x + 2/4 * y + 3/4 * z,

      p2(x) = 1/3 \
            + 1/4 * x + 2/4 * y + 3/4 * z,
            + 1/4 * x^2 + 1/2 * x*y + 3/4 * y^2 + ...

    etc.
    '''
    def f(x):
        out = 0.0
        for k in range(degree+1):
            i = 0
            tuples = _get_sum_tuples(k)
            for tup in tuples:
                # This relies on Python's 0.0**0=1.
                alpha = (i+1) / float(k+2) * x[0]**(k-i) * x[1]**i
                out += alpha
        return out

    return f


def test_generator():
    tetrahedron = numpy.array([
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 2.0],
        ])
    schemes = [
        quadrature.tetrahedron.Keast(0),
        quadrature.tetrahedron.Keast(1),
        quadrature.tetrahedron.Keast(2),
        quadrature.tetrahedron.Keast(3),
        quadrature.tetrahedron.Keast(4),
        ]
    for scheme in schemes:
        yield check_tetrahedron_scheme, scheme, tetrahedron


def check_tetrahedron_scheme(scheme, tetrahedron):
    f = _create_test_polynomial(degree=scheme.degree)
    exact_val = _integrate_exact(f, tetrahedron)
    val = quadrature.tetrahedron.integrate(f, tetrahedron, scheme)
    numpy.testing.assert_allclose(val, exact_val)
    return
