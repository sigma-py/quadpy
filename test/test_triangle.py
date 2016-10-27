# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents2

import math
import numpy
import numpy.testing
import pytest
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
    xi = sympy.DeferredVector('xi')
    x_xi = \
        + triangle[0] * (1 - xi[0] - xi[1]) \
        + triangle[1] * xi[0] \
        + triangle[2] * xi[1]
    abs_det_J = 2 * quadrature.triangle.volume(triangle)
    exact = sympy.integrate(
        sympy.integrate(abs_det_J * f(x_xi), (xi[1], 0, 1-xi[0])),
        (xi[0], 0, 1)
        )
    return float(exact)


def _integrate_monomial_over_standard_triangle(k):
    '''The integral of monomials over the standard triangle is given by

    \int_T x_0^k0 * x1^k1 = (k0!*k1!) / (2+k0+k1)!,

    see, e.g.,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Linbo Zhang, Tao Cui and Hui Liu,
    Journal of Computational Mathematics,
    Vol. 27, No. 1 (January 2009), pp. 89-96
    '''
    # exp-log to account for large values in numerator and denominator
    return math.exp(
        math.fsum([math.lgamma(kk+1) for kk in k])
        - math.lgamma(3 + sum(k))
        )


@pytest.mark.parametrize(
    'scheme',
    [quadrature.triangle.Centroid()]
    + [quadrature.triangle.Vertex()]
    + [quadrature.triangle.SevenPoint()]
    + [quadrature.triangle.HammerMarloweStroud(k) for k in range(1, 6)]
    + [quadrature.triangle.NewtonCotesClosed(k) for k in range(1, 6)]
    + [quadrature.triangle.NewtonCotesOpen(k) for k in range(6)]
    + [quadrature.triangle.Strang(k) for k in range(1, 11)]
    + [quadrature.triangle.LynessJespersen(k) for k in range(1, 22)]
    + [quadrature.triangle.Hillion(k) for k in range(1, 6)]
    + [quadrature.triangle.LaursenGellert(key) for key in [
        '1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
        '13', '14', '15a', '15b'
        ]]
    + [quadrature.triangle.Cubtri()]
    + [quadrature.triangle.Triex(19), quadrature.triangle.Triex(28)]
    + [quadrature.triangle.Dunavant(k) for k in range(1, 21)]
    + [quadrature.triangle.CoolsHaegemans(k) for k in [1, 2]]
    + [quadrature.triangle.BerntsenEspelid(k) for k in range(1, 5)]
    + [quadrature.triangle.LiuVinokur(k) for k in range(1, 14)]
    + [quadrature.triangle.WandzuraXiao(k) for k in range(1, 7)]
    + [quadrature.triangle.TaylorWingateBos(k) for k in [1, 2, 4, 5, 8]]
    + [quadrature.triangle.ZhangCuiLiu(k) for k in [1, 2, 3]]
    + [quadrature.triangle.XiaoGimbutas(k) for k in range(1, 51)]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    triangle = numpy.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0]
        ])
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in create_monomial_exponents2(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1]
            # exact_val = _integrate_exact(poly, triangle)
            exact_val = _integrate_monomial_over_standard_triangle(k)
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
    assert degree-1 >= scheme.degree
    # numpy.testing.assert_equal(degree-1, scheme.degree)
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
        # quadrature.triangle.Dunavant(20)
        # quadrature.triangle.CoolsHaegemans(2)
        quadrature.triangle.XiaoGimbutas(50)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
    scheme = quadrature.triangle.XiaoGimbutas(50)
    test_scheme(scheme)
