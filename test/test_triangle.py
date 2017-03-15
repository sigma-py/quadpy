# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents2, check_degree

import math
from matplotlib import pyplot as plt
import numpy
import pytest
import quadpy
import sympy


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
    abs_det_J = 2 * quadpy.triangle.volume(triangle)
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
    [quadpy.triangle.Centroid()]
    + [quadpy.triangle.Vertex()]
    + [quadpy.triangle.SevenPoint()]
    + [quadpy.triangle.HammerMarloweStroud(k) for k in range(1, 6)]
    + [quadpy.triangle.NewtonCotesClosed(k) for k in range(1, 6)]
    + [quadpy.triangle.NewtonCotesOpen(k) for k in range(6)]
    + [quadpy.triangle.Strang(k) for k in range(1, 11)]
    + [quadpy.triangle.LynessJespersen(k) for k in range(1, 22)]
    + [quadpy.triangle.Hillion(k) for k in range(1, 4)]
    + [quadpy.triangle.LaursenGellert(key) for key in [
        '1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
        '13', '14', '15a', '15b'
        ]]
    + [quadpy.triangle.Cubtri()]
    + [quadpy.triangle.Triex(19), quadpy.triangle.Triex(28)]
    + [quadpy.triangle.Dunavant(k) for k in range(1, 21)]
    + [quadpy.triangle.CoolsHaegemans(k) for k in [1]]
    + [quadpy.triangle.BerntsenEspelid(k) for k in range(1, 5)]
    + [quadpy.triangle.LiuVinokur(k) for k in range(1, 14)]
    + [quadpy.triangle.WandzuraXiao(k) for k in range(1, 7)]
    + [quadpy.triangle.TaylorWingateBos(k) for k in [1, 2, 4, 5, 8]]
    + [quadpy.triangle.ZhangCuiLiu(k) for k in [1, 2, 3]]
    + [quadpy.triangle.XiaoGimbutas(k) for k in range(1, 51)]
    )
def test_scheme(scheme):
    triangle = numpy.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0]
        ])
    degree = check_degree(
            lambda poly: quadpy.triangle.integrate(poly, triangle, scheme),
            _integrate_monomial_over_standard_triangle,
            create_monomial_exponents2,
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.triangle.XiaoGimbutas(10)]
    )
def test_show(scheme):
    triangle = numpy.array([
        [numpy.cos(0.5*numpy.pi), numpy.sin(0.5*numpy.pi)],
        [numpy.cos(7.0/6.0*numpy.pi), numpy.sin(7.0/6.0*numpy.pi)],
        [numpy.cos(11.0/6.0*numpy.pi), numpy.sin(11.0/6.0*numpy.pi)],
        ])
    quadpy.triangle.show(scheme, triangle)
    return


if __name__ == '__main__':
    scheme = quadpy.triangle.XiaoGimbutas(50)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
