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


def _integrate_exact(f, triangle):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference triangle [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the triangle. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    def g(xi):
        pxi = triangle[0] * (1 - xi[0] - xi[1]) \
            + triangle[1] * xi[0] \
            + triangle[2] * xi[1]
        return f(pxi)

    x = sympy.DeferredVector('x')
    exact = 2 * quadrature.triangle.volume(triangle) \
        * sympy.integrate(
            sympy.integrate(g(x), (x[1], 0, 1-x[0])),
            (x[0], 0, 1)
            )
    return float(exact)


def _create_test_polynomial(degree):
    '''The k-th order terms in polynomial have the form

        alpha_{k,i} * x^{k-i} * y^k,

    with i in {0,...,k}. Take

      alpha_{k, i} = (i+1) / (k+2)

    such that

      p0(x) = 1/2,

      p1(x) = 1/2 \
            + 1/3 * x + 2/3 * y,

      p1(x) = 1/2 \
            + 1/3 * x + 2/3 * y \
            + 1/4 * x^2 + 1/2 * x*y + 3/4 * y^2

    etc.
    '''
    def f(x):
        out = 0.0
        for k in range(degree+1):
            i = 0
            for i in range(k+1):
                # This relies on Python's 0.0**0=1.
                alpha = (i+1) / float(k+2) * x[0]**(k-i) * x[1]**i
                out += alpha
        return out

    return f


def test_generator():
    triangle = numpy.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.6, 0.5]
        ])
    schemes = [
        quadrature.triangle.Centroid(),
        quadrature.triangle.Vertex(),
        quadrature.triangle.SevenPoint(),
        quadrature.triangle.Strang1(),
        quadrature.triangle.Strang2(),
        quadrature.triangle.Strang3(),
        quadrature.triangle.Strang4(),
        quadrature.triangle.Strang5(),
        quadrature.triangle.Strang6(),
        quadrature.triangle.Strang7(),
        quadrature.triangle.Strang8(),
        quadrature.triangle.Strang9(),
        quadrature.triangle.Strang10(),
        quadrature.triangle.Toms584_19(),
        quadrature.triangle.Toms612_19(),
        quadrature.triangle.Toms612_28(),
        quadrature.triangle.Toms706_37(),
        quadrature.triangle.Gauss4x4(),
        quadrature.triangle.Gauss8x8(),
        quadrature.triangle.Dunavant(1),
        quadrature.triangle.Dunavant(2),
        quadrature.triangle.Dunavant(3),
        quadrature.triangle.Dunavant(4),
        quadrature.triangle.Dunavant(5),
        quadrature.triangle.Dunavant(6),
        quadrature.triangle.Dunavant(7),
        quadrature.triangle.Dunavant(8),
        quadrature.triangle.Dunavant(9),
        quadrature.triangle.Dunavant(10),
        quadrature.triangle.Dunavant(11),
        quadrature.triangle.Dunavant(12),
        quadrature.triangle.Dunavant(13),
        quadrature.triangle.Dunavant(14),
        quadrature.triangle.Dunavant(15),
        ]
    for scheme in schemes:
        yield check_triangle_scheme, scheme, triangle


def check_triangle_scheme(scheme, triangle):
    f = _create_test_polynomial(degree=scheme.degree)
    exact_val = _integrate_exact(f, triangle)
    val = quadrature.triangle.integrate(f, triangle, scheme)
    numpy.testing.assert_allclose(val, exact_val)
    return


def test_show():
    triangle = numpy.array([
        [numpy.cos(0.5*numpy.pi), numpy.sin(0.5*numpy.pi)],
        [numpy.cos(7.0/6.0*numpy.pi), numpy.sin(7.0/6.0*numpy.pi)],
        [numpy.cos(11.0/6.0*numpy.pi), numpy.sin(11.0/6.0*numpy.pi)],
        ])
    quadrature.triangle.show(
        triangle,
        # quadrature.triangle.Centroid()
        # quadrature.triangle.Vertex()
        # quadrature.triangle.SevenPoint()
        # quadrature.triangle.Strang9()
        quadrature.triangle.Dunavant(12)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
