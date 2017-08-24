# -*- coding: utf-8 -*-
#
from helpers import (
    check_degree, integrate_monomial_over_standard_simplex
    )

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
    'scheme,tol',
    [(quadpy.triangle.BerntsenEspelid(k), 1.0e-14) for k in range(1, 5)]
    + [(quadpy.triangle.Centroid(), 1.0e-14)]
    + [(quadpy.triangle.CoolsHaegemans(k), 1.0e-14) for k in [1]]
    + [(quadpy.triangle.Cubtri(), 1.0e-14)]
    + [(quadpy.triangle.Dunavant(k), 1.0e-14) for k in range(1, 21)]
    + [(quadpy.triangle.Gatermann(), 1.0e-14)]
    + [(quadpy.triangle.GrundmannMoeller(k), 1.0e-14) for k in range(10)]
    + [(quadpy.triangle.HammerMarloweStroud(k), 1.0e-14) for k in range(1, 6)]
    + [(quadpy.triangle.HammerStroud(k), 1.0e-14) for k in [2, 3]]
    + [(quadpy.triangle.Hillion(k), 1.0e-14) for k in range(1, 4)]
    + [(quadpy.triangle.Hillion(k), 1.0e-6) for k in [4, 5, 6, 8, 10]]
    + [(quadpy.triangle.LaursenGellert(key), 1.0e-14) for key in [
        '1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
        '13', '14', '15a', '15b'
        ]]
    + [(quadpy.triangle.Lether(k), 1.0e-14) for k in range(1, 14)]
    + [(quadpy.triangle.LiuVinokur(k), 1.0e-14) for k in range(1, 14)]
    + [(quadpy.triangle.LynessJespersen(k), 1.0e-14) for k in range(1, 22)]
    + [(quadpy.triangle.NewtonCotesClosed(k), 1.0e-14) for k in range(1, 6)]
    + [(quadpy.triangle.NewtonCotesOpen(k), 1.0e-14) for k in range(6)]
    + [(quadpy.triangle.Papanicolopulos('fs', k), 1.0e-14) for k in range(9)]
    + [(quadpy.triangle.Papanicolopulos('rot', k), 1.0e-14) for k in range(18)]
    + [(quadpy.triangle.SevenPoint(), 1.0e-14)]
    + [(quadpy.triangle.Strang(k), 1.0e-14) for k in range(1, 11)]
    + [(quadpy.triangle.Stroud(k), 1.0e-14) for k in [
        'T2 3-1', 'T2 5-1', 'T2 7-1'
        ]
       ]
    + [(quadpy.triangle.TaylorWingateBos(k), 1.0e-14) for k in [1, 2, 4, 5, 8]]
    + [(quadpy.triangle.Triex(k), 1.0e-14) for k in [19, 28]]
    + [(quadpy.triangle.Vertex(), 1.0e-14)]
    + [(quadpy.triangle.VioreanuRokhlin(k), 1.0e-14) for k in range(20)]
    + [(quadpy.triangle.Walkington(k), 1.0e-14) for k in [1, 2, 3, 5, 'p5']]
    + [(quadpy.triangle.WandzuraXiao(k), 1.0e-14) for k in range(1, 7)]
    + [(quadpy.triangle.WilliamsShunnJameson(k), 1.0e-14) for k in range(1, 9)]
    + [(quadpy.triangle.WitherdenVincent(k), 1.0e-14) for k in [
        1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
        ]]
    + [(quadpy.triangle.XiaoGimbutas(k), 1.0e-14) for k in range(1, 51)]
    + [(quadpy.triangle.ZhangCuiLiu(k), 1.0e-14) for k in [1, 2, 3]]
    )
def test_scheme(scheme, tol):
    triangle = numpy.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0]
        ])
    degree = check_degree(
            lambda poly: quadpy.triangle.integrate(poly, triangle, scheme),
            integrate_monomial_over_standard_simplex,
            lambda n: quadpy.helpers.partition(n, 2),
            scheme.degree + 1,
            tol=tol
            )
    assert degree >= scheme.degree, \
        'Observed: {}, expected: {}'.format(degree, scheme.degree)
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
    scheme_ = quadpy.triangle.Stroud('T2 7-1')
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
