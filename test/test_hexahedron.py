# -*- coding: utf-8 -*-
#
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy
import pytest
import sympy

import quadpy

from helpers import check_degree

schemes = (
    [quadpy.hexahedron.product(quadpy.line_segment.Midpoint())]
    + [quadpy.hexahedron.product(quadpy.line_segment.Trapezoidal())]
    + [
        quadpy.hexahedron.product(quadpy.line_segment.GaussLegendre(k))
        for k in range(1, 6)
    ]
    + [
        quadpy.hexahedron.product(quadpy.line_segment.NewtonCotesClosed(k))
        for k in range(1, 5)
    ]
    + [
        quadpy.hexahedron.product(quadpy.line_segment.NewtonCotesOpen(k))
        for k in range(5)
    ]
    + [
        quadpy.hexahedron.hammer_stroud_1_3(),
        quadpy.hexahedron.hammer_stroud_2_3(),
        quadpy.hexahedron.hammer_stroud_4_3(),
        quadpy.hexahedron.hammer_stroud_5_3a(),
        quadpy.hexahedron.hammer_stroud_5_3b(),
        quadpy.hexahedron.hammer_stroud_6_3(),
        quadpy.hexahedron.hammer_wymore(),
        quadpy.hexahedron.mustard_lyness_blatt_1(),
        quadpy.hexahedron.mustard_lyness_blatt_2(),
        quadpy.hexahedron.mustard_lyness_blatt_3(),
        quadpy.hexahedron.mustard_lyness_blatt_4(),
        quadpy.hexahedron.mustard_lyness_blatt_5(),
        quadpy.hexahedron.mustard_lyness_blatt_6(),
        quadpy.hexahedron.mustard_lyness_blatt_7(),
        quadpy.hexahedron.sadowsky(),
        quadpy.hexahedron.stroud_c3_3_1(),
        quadpy.hexahedron.stroud_c3_3_2(),
        quadpy.hexahedron.stroud_c3_3_3(),
        quadpy.hexahedron.stroud_c3_3_4(),
        quadpy.hexahedron.stroud_c3_3_5(),
        quadpy.hexahedron.stroud_c3_3_6(),
        quadpy.hexahedron.stroud_c3_3_7(),
        quadpy.hexahedron.stroud_c3_5_1(),
        quadpy.hexahedron.stroud_c3_5_2(),
        quadpy.hexahedron.stroud_c3_5_3(),
        quadpy.hexahedron.stroud_c3_5_4(),
        quadpy.hexahedron.stroud_c3_5_5(),
        quadpy.hexahedron.stroud_c3_5_6(),
        quadpy.hexahedron.stroud_c3_5_7(),
        quadpy.hexahedron.stroud_c3_5_8(),
        quadpy.hexahedron.stroud_c3_7_1a(),
        quadpy.hexahedron.stroud_c3_7_1b(),
        quadpy.hexahedron.stroud_c3_7_2(),
        quadpy.hexahedron.stroud_c3_7_3(),
        quadpy.hexahedron.stroud_1967(),
        quadpy.hexahedron.tyler_1(),
        quadpy.hexahedron.tyler_2(),
    ]
    + [
        # quadpy.ncube.dobrodeev_1970(3),
        quadpy.ncube.dobrodeev_1978(3),
        quadpy.ncube.hammer_stroud_1n(3),
        quadpy.ncube.hammer_stroud_2n(3),
        quadpy.ncube.stroud_cn_1_1(3),
        quadpy.ncube.stroud_cn_1_2(3),
        quadpy.ncube.stroud_cn_2_1(3),
        quadpy.ncube.stroud_cn_2_2(3),
        quadpy.ncube.stroud_cn_3_1(3),
        quadpy.ncube.stroud_cn_3_2(3),
        quadpy.ncube.stroud_cn_3_3(3),
        quadpy.ncube.stroud_cn_3_4(3),
        quadpy.ncube.stroud_cn_3_5(3),
        quadpy.ncube.stroud_cn_3_6(3),
        quadpy.ncube.stroud_cn_5_2(3),
        quadpy.ncube.stroud_cn_5_3(3),
        quadpy.ncube.stroud_cn_5_4(3),
        quadpy.ncube.stroud_cn_5_5(3),
        quadpy.ncube.stroud_cn_5_6(3),
        quadpy.ncube.stroud_cn_5_7(3),
        quadpy.ncube.stroud_cn_5_8(3),
        quadpy.ncube.stroud_cn_5_9(3),
        quadpy.ncube.stroud_cn_7_1(3),
    ]
)


def _integrate_exact(f, hexa):
    xi = sympy.DeferredVector("xi")
    pxi = (
        +hexa[0] * 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
        + hexa[1] * 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
        + hexa[2] * 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
        + hexa[3] * 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
        + hexa[4] * 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
        + hexa[5] * 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
        + hexa[6] * 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])
        + hexa[7] * 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])
    )
    pxi = [sympy.expand(pxi[0]), sympy.expand(pxi[1]), sympy.expand(pxi[2])]
    # determinant of the transformation matrix
    J = sympy.Matrix(
        [
            [
                sympy.diff(pxi[0], xi[0]),
                sympy.diff(pxi[0], xi[1]),
                sympy.diff(pxi[0], xi[2]),
            ],
            [
                sympy.diff(pxi[1], xi[0]),
                sympy.diff(pxi[1], xi[1]),
                sympy.diff(pxi[1], xi[2]),
            ],
            [
                sympy.diff(pxi[2], xi[0]),
                sympy.diff(pxi[2], xi[1]),
                sympy.diff(pxi[2], xi[2]),
            ],
        ]
    )
    det_J = sympy.det(J)
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.

    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))
    g_xi = f(pxi)
    exact = sympy.integrate(
        sympy.integrate(
            sympy.integrate(abs_det_J * g_xi, (xi[2], -1, 1)), (xi[1], -1, 1)
        ),
        (xi[0], -1, 1),
    )
    return float(exact)


def _integrate_exact2(k, x0, x1, y0, y1, z0, z1):
    return (
        1.0
        / (k[0] + 1)
        * (x1 ** (k[0] + 1) - x0 ** (k[0] + 1))
        * 1.0
        / (k[1] + 1)
        * (y1 ** (k[1] + 1) - y0 ** (k[1] + 1))
        * 1.0
        / (k[2] + 1)
        * (z1 ** (k[2] + 1) - z0 ** (k[2] + 1))
    )


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme, tol=1.0e-14, print_degree=False):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    x = [-1.0, +1.0]
    y = [-1.0, +1.0]
    z = [-1.0, +1.0]
    hexa = quadpy.hexahedron.cube_points(x, y, z)
    degree = check_degree(
        lambda poly: scheme.integrate(poly, hexa),
        lambda k: _integrate_exact2(k, x[0], x[1], y[0], y[1], z[0], z[1]),
        3,
        scheme.degree + 1,
        tol=tol,
    )
    if print_degree:
        print("Detected degree {}, scheme degree {}.".format(degree, scheme.degree))
    assert degree == scheme.degree, scheme.name
    return


@pytest.mark.parametrize(
    "scheme", [quadpy.hexahedron.product(quadpy.line_segment.NewtonCotesClosed(2))]
)
def test_show(scheme):
    scheme.show(backend="mpl")
    plt.close()
    return


if __name__ == "__main__":
    # scheme_ = Product(quadpy.line_segment.NewtonCotesOpen(5))
    # scheme_ = quadpy.hexahedron.HammerStroud("6-3")
    # test_scheme(scheme_, 1.0e-14, print_degree=True)
    # test_show(scheme_)
    # scheme_.show(backend="vtk")
    from helpers import find_equal

    find_equal(schemes)
