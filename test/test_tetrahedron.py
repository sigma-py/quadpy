# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy
import pytest
import sympy

import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex

from helpers import check_degree


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
    xi = sympy.DeferredVector("xi")
    x_xi = (
        +tetrahedron[0] * (1 - xi[0] - xi[1] - xi[2])
        + tetrahedron[1] * xi[0]
        + tetrahedron[2] * xi[1]
        + tetrahedron[3] * xi[2]
    )
    abs_det_J = 6 * quadpy.tetrahedron.volume(tetrahedron)
    exact = sympy.integrate(
        sympy.integrate(
            sympy.integrate(abs_det_J * f(x_xi), (xi[2], 0, 1 - xi[0] - xi[1])),
            (xi[1], 0, 1 - xi[0]),
        ),
        (xi[0], 0, 1),
    )
    return float(exact)


@pytest.mark.parametrize(
    "scheme",
    [quadpy.tetrahedron.BeckersHaegemans(k) for k in [8, 9]]
    + [quadpy.tetrahedron.Gatermann()]
    + [quadpy.tetrahedron.GrundmannMoeller(k) for k in range(8)]
    + [quadpy.tetrahedron.HammerStroud(k) for k in [2, 3]]
    + [quadpy.tetrahedron.HammerMarloweStroud(k) for k in [1, 2, 3]]
    + [quadpy.tetrahedron.Keast(k) for k in range(11)]
    + [quadpy.tetrahedron.LiuVinokur(k) for k in range(1, 15)]
    + [quadpy.tetrahedron.MaeztuSainz()]
    + [quadpy.tetrahedron.NewtonCotesClosed(k) for k in range(1, 7)]
    + [quadpy.tetrahedron.NewtonCotesOpen(k) for k in range(7)]
    + [quadpy.tetrahedron.ShunnHam(k) for k in range(1, 7)]
    + [quadpy.tetrahedron.Stroud(k) for k in ["T3 5-1", "T3 7-1"]]
    + [quadpy.tetrahedron.VioreanuRokhlin(k) for k in range(10)]
    + [quadpy.tetrahedron.Walkington(k) for k in [1, 2, 3, 5, "p5", 7]]
    + [quadpy.tetrahedron.WilliamsShunnJameson()]
    + [quadpy.tetrahedron.WitherdenVincent(k) for k in [1, 2, 3, 5, 6, 7, 8, 9, 10]]
    + [quadpy.tetrahedron.XiaoGimbutas(k) for k in range(1, 16)]
    + [quadpy.tetrahedron.Yu(k) for k in range(1, 6)]
    + [quadpy.tetrahedron.ZhangCuiLiu(k) for k in [1, 2]]
    + [quadpy.tetrahedron.Zienkiewicz(k) for k in [4, 5]],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    tetrahedron = numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    degree = check_degree(
        lambda poly: quadpy.tetrahedron.integrate(poly, tetrahedron, scheme),
        integrate_monomial_over_unit_simplex,
        3,
        scheme.degree + 1,
    )
    assert degree == scheme.degree, "Observed: {}, expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.tetrahedron.HammerMarloweStroud(3)])
def test_show(scheme):
    tet = numpy.array(
        [
            [numpy.cos(0.5 * numpy.pi), numpy.sin(0.5 * numpy.pi), -0.5],
            [numpy.cos(7.0 / 6.0 * numpy.pi), numpy.sin(7.0 / 6.0 * numpy.pi), -0.5],
            [numpy.cos(11.0 / 6.0 * numpy.pi), numpy.sin(11.0 / 6.0 * numpy.pi), -0.5],
            [0.0, 0.0, 1.0],
        ]
    )
    quadpy.tetrahedron.show(scheme, tet)
    plt.close()
    return


if __name__ == "__main__":
    scheme_ = quadpy.tetrahedron.Stroud("T3 7-1")
    test_scheme(scheme_)
    # test_show(scheme_)
    quadpy.tetrahedron.show(scheme_, backend="vtk")
