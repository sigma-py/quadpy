# -*- coding: utf-8 -*-
#
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy
import pytest
import sympy

import quadpy
from quadpy.hexahedron import Product

from helpers import check_degree


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


@pytest.mark.parametrize(
    "scheme, tol",
    [(Product(quadpy.line_segment.Midpoint()), 1.0e-14)]
    + [(Product(quadpy.line_segment.Trapezoidal()), 1.0e-14)]
    + [(Product(quadpy.line_segment.GaussLegendre(k)), 1.0e-14) for k in range(1, 6)]
    + [
        (Product(quadpy.line_segment.NewtonCotesClosed(k)), 1.0e-14)
        for k in range(1, 5)
    ]
    + [(Product(quadpy.line_segment.NewtonCotesOpen(k)), 1.0e-14) for k in range(5)]
    + [
        (quadpy.hexahedron.HammerStroud[k](), 1.0e-14)
        for k in ["1-3", "2-3", "4-3", "5-3a", "5-3b", "6-3"]
    ]
    + [
        (quadpy.hexahedron.Stroud[k](), 1.0e-14)
        for k in [
            "C3 3-1",
            "C3 3-2",
            "C3 3-3",
            "C3 3-4",
            "C3 3-5",
            "C3 3-6",
            "C3 3-7",
            "C3 5-1",
            "C3 5-2",
            "C3 5-3",
            "C3 5-4",
            "C3 5-5",
            "C3 5-6",
            "C3 5-7",
            "C3 5-8",
            "C3 7-1a",
            "C3 7-1b",
            "C3 7-2",
            "C3 7-3",
        ]
    ]
    + [
        (quadpy.hexahedron.StroudN(k), 1.0e-14)
        for k in [
            "Cn 1-1",
            "Cn 1-2",
            "Cn 2-1",
            "Cn 2-2",
            "Cn 3-1",
            "Cn 3-2",
            "Cn 3-3",
            "Cn 3-4",
            "Cn 3-5",
            "Cn 3-6",
            "Cn 5-2",
            "Cn 5-3",
            "Cn 5-4",
            "Cn 5-5",
            "Cn 5-6",
            "Cn 5-7",
            "Cn 5-8",
            "Cn 5-9",
        ]
    ]
    + [(quadpy.hexahedron.StroudN(k), 1.0e-7) for k in ["Cn 7-1"]],
)
def test_scheme(scheme, tol, print_degree=False):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    x = [-1.0, +1.0]
    y = [-1.0, +1.0]
    z = [-1.0, +1.0]
    hexa = quadpy.hexahedron.cube_points(x, y, z)
    degree = check_degree(
        lambda poly: quadpy.hexahedron.integrate(poly, hexa, scheme),
        lambda k: _integrate_exact2(k, x[0], x[1], y[0], y[1], z[0], z[1]),
        3,
        scheme.degree + 1,
        tol=tol,
    )
    if print_degree:
        print("Detected degree {}, scheme degree {}.".format(degree, scheme.degree))
    assert degree == scheme.degree
    return


@pytest.mark.parametrize("scheme", [Product(quadpy.line_segment.NewtonCotesClosed(2))])
def test_show(scheme):
    quadpy.hexahedron.show(scheme)
    plt.close()
    return


if __name__ == "__main__":
    # scheme_ = Product(quadpy.line_segment.NewtonCotesOpen(5))
    scheme_ = quadpy.hexahedron.HammerStroud("6-3")
    test_scheme(scheme_, 1.0e-14, print_degree=True)
    # test_show(scheme_)
    quadpy.hexahedron.show(scheme_, backend="vtk")
