# -*- coding: utf-8 -*-
#
from helpers import (
    partition, check_degree, integrate_monomial_over_standard_simplex
    )

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


@pytest.mark.parametrize(
    'scheme',
    [quadpy.triangle.BerntsenEspelid(k) for k in range(1, 5)]
    + [quadpy.triangle.Centroid()]
    + [quadpy.triangle.CoolsHaegemans(k) for k in [1]]
    + [quadpy.triangle.Cubtri()]
    + [quadpy.triangle.Dunavant(k) for k in range(1, 21)]
    + [quadpy.triangle.Gatermann()]
    + [quadpy.triangle.GrundmannMoeller(k) for k in range(10)]
    + [quadpy.triangle.HammerMarloweStroud(k) for k in range(1, 6)]
    + [quadpy.triangle.Hillion(k) for k in range(1, 4)]
    + [quadpy.triangle.LaursenGellert(key) for key in [
        '1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
        '13', '14', '15a', '15b'
        ]]
    + [quadpy.triangle.LiuVinokur(k) for k in range(1, 14)]
    + [quadpy.triangle.LynessJespersen(k) for k in range(1, 22)]
    + [quadpy.triangle.NewtonCotesClosed(k) for k in range(1, 6)]
    + [quadpy.triangle.NewtonCotesOpen(k) for k in range(6)]
    + [quadpy.triangle.SevenPoint()]
    + [quadpy.triangle.Strang(k) for k in range(1, 11)]
    + [quadpy.triangle.Stroud(k) for k in range(10)]
    + [quadpy.triangle.TaylorWingateBos(k) for k in [1, 2, 4, 5, 8]]
    + [quadpy.triangle.Triex(19), quadpy.triangle.Triex(28)]
    + [quadpy.triangle.Vertex()]
    + [quadpy.triangle.VioreanuRokhlin(k) for k in range(20)]
    + [quadpy.triangle.Walkington(k) for k in [1, 2, 3, 5, 'p5']]
    + [quadpy.triangle.WandzuraXiao(k) for k in range(1, 7)]
    + [quadpy.triangle.WilliamsShunnJameson(k) for k in range(1, 9)]
    + [quadpy.triangle.XiaoGimbutas(k) for k in range(1, 51)]
    + [quadpy.triangle.ZhangCuiLiu(k) for k in [1, 2, 3]]
    )
def test_scheme(scheme):
    triangle = numpy.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0]
        ])
    degree = check_degree(
            lambda poly: quadpy.triangle.integrate(poly, triangle, scheme),
            integrate_monomial_over_standard_simplex,
            lambda n: partition(n, 2),
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


def test_volume():
    # Assert computation of triangle volume in 3D is correct
    triangle = numpy.array([
        [0.0, 0.0, 0.0],
        [1.0, 2.0, 3.0],
        [0.7, 0.4, 1.1],
        ])
    ref = numpy.sqrt(3.0) / 2.0
    assert abs(quadpy.triangle.get_vol(triangle) - ref) < 1.0e-14 * ref

    triangle = numpy.array([
        [0.0, 0.0, 0.0],
        [0.3, 0.4, 0.5],
        [0.7, 0.4, 1.1],
        ])
    ref = numpy.sqrt(0.0209)
    assert abs(quadpy.triangle.get_vol(triangle) - ref) < 1.0e-14 * ref
    return


if __name__ == '__main__':
    scheme_ = quadpy.triangle.VioreanuRokhlin(10)
    test_scheme(scheme_)
    test_show(scheme_)
    plt.show()
