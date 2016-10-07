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


def _create_test_monomials(degree):
    '''Returns a list of all monomials of degree :degree:.
    '''
    return [
        lambda x: x[0]**(degree-k) * x[1]**k
        for k in range(degree+1)
        ]


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
        quadrature.triangle.Strang(1),
        quadrature.triangle.Strang(2),
        quadrature.triangle.Strang(3),
        quadrature.triangle.Strang(4),
        quadrature.triangle.Strang(5),
        quadrature.triangle.Strang(6),
        quadrature.triangle.Strang(7),
        quadrature.triangle.Strang(8),
        quadrature.triangle.Strang(9),
        quadrature.triangle.Strang(10),
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
        quadrature.triangle.Dunavant(16),
        quadrature.triangle.Dunavant(17),
        # quadrature.triangle.Dunavant(18),
        quadrature.triangle.Dunavant(19),
        quadrature.triangle.Dunavant(20),
        ]
    for scheme in schemes:
        yield check_scheme, scheme, triangle


def check_scheme(scheme, triangle):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    success = True
    degree = 0
    max_degree = 100
    while success:
        for poly in _create_test_monomials(degree):
            exact_val = _integrate_exact(poly, triangle)
            val = quadrature.triangle.integrate(
                    poly, triangle, scheme
                    )
            if abs(exact_val - val) > 1.0e-10:
                success = False
                break
        if not success:
            break
        if degree >= max_degree:
            break
        degree += 1
    numpy.testing.assert_equal(degree-1, scheme.degree)
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
        # quadrature.triangle.Strang(9)
        quadrature.triangle.Dunavant(20)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
