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


@pytest.mark.parametrize('scheme', [
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
    quadrature.triangle.Dunavant(18),
    quadrature.triangle.Dunavant(19),
    quadrature.triangle.Dunavant(20),
    quadrature.triangle.ZhangCuiLiu(1),
    quadrature.triangle.ZhangCuiLiu(2),
    quadrature.triangle.ZhangCuiLiu(3),
    quadrature.triangle.WandzuraXiao(1),
    quadrature.triangle.WandzuraXiao(2),
    quadrature.triangle.WandzuraXiao(3),
    quadrature.triangle.WandzuraXiao(4),
    quadrature.triangle.WandzuraXiao(5),
    quadrature.triangle.WandzuraXiao(6),
    quadrature.triangle.LynessJespersen(1),
    quadrature.triangle.LynessJespersen(2),
    quadrature.triangle.LynessJespersen(3),
    quadrature.triangle.LynessJespersen(4),
    quadrature.triangle.LynessJespersen(5),
    quadrature.triangle.LynessJespersen(6),
    quadrature.triangle.LynessJespersen(7),
    quadrature.triangle.LynessJespersen(8),
    quadrature.triangle.LynessJespersen(9),
    quadrature.triangle.LynessJespersen(10),
    quadrature.triangle.LynessJespersen(11),
    quadrature.triangle.LynessJespersen(12),
    quadrature.triangle.LynessJespersen(13),
    quadrature.triangle.LynessJespersen(14),
    quadrature.triangle.LynessJespersen(15),
    quadrature.triangle.LynessJespersen(16),
    quadrature.triangle.LynessJespersen(17),
    quadrature.triangle.LynessJespersen(18),
    quadrature.triangle.LynessJespersen(19),
    quadrature.triangle.LynessJespersen(20),
    quadrature.triangle.LynessJespersen(21),
    quadrature.triangle.NewtonCotesClosed(1),
    quadrature.triangle.NewtonCotesClosed(2),
    quadrature.triangle.NewtonCotesClosed(3),
    quadrature.triangle.NewtonCotesClosed(4),
    quadrature.triangle.NewtonCotesClosed(5),
    quadrature.triangle.NewtonCotesOpen(0),
    quadrature.triangle.NewtonCotesOpen(1),
    quadrature.triangle.NewtonCotesOpen(2),
    quadrature.triangle.NewtonCotesOpen(3),
    quadrature.triangle.NewtonCotesOpen(4),
    quadrature.triangle.NewtonCotesOpen(5),
    quadrature.triangle.TaylorWingateBos(1),
    quadrature.triangle.TaylorWingateBos(2),
    quadrature.triangle.TaylorWingateBos(4),
    quadrature.triangle.TaylorWingateBos(5),
    quadrature.triangle.TaylorWingateBos(8),
    quadrature.triangle.BerntsenEspelid(1),
    quadrature.triangle.BerntsenEspelid(2),
    quadrature.triangle.BerntsenEspelid(3),
    quadrature.triangle.BerntsenEspelid(4),
    quadrature.triangle.HammerMarloweStroud(1),
    quadrature.triangle.HammerMarloweStroud(2),
    quadrature.triangle.HammerMarloweStroud(3),
    quadrature.triangle.HammerMarloweStroud(4),
    quadrature.triangle.HammerMarloweStroud(5),
    quadrature.triangle.Cowper(1),
    quadrature.triangle.Cowper(2),
    quadrature.triangle.Cowper(3),
    quadrature.triangle.Cowper(4),
    quadrature.triangle.Cowper(5),
    quadrature.triangle.Cowper(6),
    quadrature.triangle.Cowper(7),
    quadrature.triangle.Cowper(8),
    quadrature.triangle.Cowper(9),
    quadrature.triangle.Cowper(10),
    quadrature.triangle.LiuVinokur(1),
    quadrature.triangle.LiuVinokur(2),
    quadrature.triangle.LiuVinokur(3),
    quadrature.triangle.LiuVinokur(4),
    quadrature.triangle.LiuVinokur(5),
    quadrature.triangle.LiuVinokur(6),
    quadrature.triangle.LiuVinokur(7),
    quadrature.triangle.LiuVinokur(8),
    quadrature.triangle.LiuVinokur(9),
    quadrature.triangle.LiuVinokur(10),
    quadrature.triangle.LiuVinokur(11),
    quadrature.triangle.LiuVinokur(12),
    quadrature.triangle.LiuVinokur(13),
    quadrature.triangle.Hillion(1),
    quadrature.triangle.Hillion(2),
    quadrature.triangle.Hillion(3),
    quadrature.triangle.Hillion(4),
    quadrature.triangle.Hillion(5),
    quadrature.triangle.CoolsHaegemans(1),
    quadrature.triangle.CoolsHaegemans(2),
    quadrature.triangle.LaursenGellert('1'),
    quadrature.triangle.LaursenGellert('2a'),
    quadrature.triangle.LaursenGellert('2b'),
    quadrature.triangle.LaursenGellert('3'),
    quadrature.triangle.LaursenGellert('4'),
    quadrature.triangle.LaursenGellert('5'),
    quadrature.triangle.LaursenGellert('6'),
    quadrature.triangle.LaursenGellert('7'),
    quadrature.triangle.LaursenGellert('8'),
    quadrature.triangle.LaursenGellert('9'),
    quadrature.triangle.LaursenGellert('10'),
    quadrature.triangle.LaursenGellert('11'),
    quadrature.triangle.LaursenGellert('12'),
    quadrature.triangle.LaursenGellert('13'),
    quadrature.triangle.LaursenGellert('14'),
    quadrature.triangle.LaursenGellert('15a'),
    quadrature.triangle.LaursenGellert('15b'),
    ])
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
        quadrature.triangle.Dunavant(20)
        # quadrature.triangle.CoolsHaegemans(2)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
    scheme = quadrature.triangle.CoolsHaegemans(2)
    test_scheme(scheme)
