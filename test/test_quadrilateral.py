# -*- coding: utf-8 -*-
#
from __future__ import print_function

import numpy
import pytest
import sympy

import orthopy
import quadpy
from quadpy.quadrilateral import Product

from helpers import check_degree_ortho


def _integrate_exact(f, quadrilateral):
    xi = sympy.DeferredVector("xi")
    pxi = (
        quadrilateral[0] * 0.25 * (1.0 + xi[0]) * (1.0 + xi[1])
        + quadrilateral[1] * 0.25 * (1.0 - xi[0]) * (1.0 + xi[1])
        + quadrilateral[2] * 0.25 * (1.0 - xi[0]) * (1.0 - xi[1])
        + quadrilateral[3] * 0.25 * (1.0 + xi[0]) * (1.0 - xi[1])
    )
    pxi = [sympy.expand(pxi[0]), sympy.expand(pxi[1])]
    # determinant of the transformation matrix
    det_J = +sympy.diff(pxi[0], xi[0]) * sympy.diff(pxi[1], xi[1]) - sympy.diff(
        pxi[1], xi[0]
    ) * sympy.diff(pxi[0], xi[1])
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.
    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))

    g_xi = f(pxi)

    exact = sympy.integrate(
        sympy.integrate(abs_det_J * g_xi, (xi[1], -1, 1)), (xi[0], -1, 1)
    )
    return float(exact)


def _integrate_exact2(k, x0, x1, y0, y1):
    return (
        1.0
        / (k[0] + 1)
        * (x1 ** (k[0] + 1) - x0 ** (k[0] + 1))
        * 1.0
        / (k[1] + 1)
        * (y1 ** (k[1] + 1) - y0 ** (k[1] + 1))
    )


@pytest.mark.parametrize(
    "scheme,tol",
    [(quadpy.quadrilateral.CoolsHaegemans1985(k), 1.0e-10) for k in range(1, 4)]
    + [(quadpy.quadrilateral.CoolsHaegemans1988(k), 1.0e-14) for k in [1, 2]]
    + [(quadpy.quadrilateral.Dunavant(k), 1.0e-13) for k in range(11)]
    + [(quadpy.quadrilateral.HammerStroud(k), 1.0e-14) for k in ["1-2", "2-2", "3-2"]]
    + [(quadpy.quadrilateral.MorrowPatterson(k), 1.0e-5) for k in [1, 2]]
    + [
        (quadpy.quadrilateral.Stroud(k), 1.0e-13)
        for k in [
            "C2 1-1",
            "C2 1-2",
            "C2 3-1",
            "C2 3-2",
            "C2 3-3",
            "C2 3-4",
            "C2 3-5",
            "C2 5-1",
            "C2 5-2",
            "C2 5-3",
            "C2 5-4",
            "C2 5-5",
            "C2 5-6",
            "C2 5-7",
            "C2 7-1",
            "C2 7-2",
            "C2 7-3",
            "C2 7-4",
            "C2 7-5",
            "C2 7-6",
            "C2 9-1",
            "C2 11-1",
            "C2 11-2",
            "C2 13-1",
            "C2 15-1",
            "C2 15-2",
        ]
    ]
    + [
        (quadpy.quadrilateral.StroudN(k), 1.0e-14)
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
            "Cn 5-9",
        ]
    ]
    + [(quadpy.quadrilateral.HaegemansPiessens(), 1.0e-14)]
    + [(quadpy.quadrilateral.PiessensHaegemans(k), 1.0e-14) for k in [1, 2]]
    # TODO better-quality points/weights for Schmidt
    + [(quadpy.quadrilateral.Schmid(k), 1.0e-10) for k in [2, 4, 6]]
    + [(quadpy.quadrilateral.Sommariva(k), 1.0e-13) for k in range(1, 56)]
    + [(quadpy.quadrilateral.StroudN(k), 1.0e-8) for k in ["Cn 7-1"]]
    + [
        (quadpy.quadrilateral.WissmannBecker(k), 1.0e-14)
        for k in ["4-1", "4-2", "6-1", "6-2", "8-1", "8-2"]
    ]
    + [
        (quadpy.quadrilateral.WitherdenVincent(k), 1.0e-14)
        for k in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    ]
    + [(Product(quadpy.line_segment.Midpoint()), 1.0e-14)]
    + [(Product(quadpy.line_segment.Trapezoidal()), 1.0e-14)]
    + [(Product(quadpy.line_segment.GaussLegendre(k)), 1.0e-14) for k in range(1, 5)]
    + [
        (Product(quadpy.line_segment.NewtonCotesClosed(k)), 1.0e-14)
        for k in range(1, 5)
    ]
    + [(Product(quadpy.line_segment.NewtonCotesOpen(k)), 1.0e-14) for k in range(6)],
)
def test_scheme(scheme, tol):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    def eval_orthopolys(x):
        return numpy.concatenate(
            orthopy.quadrilateral.tree(x, scheme.degree + 1, symbolic=False)
        )

    quad = quadpy.quadrilateral.rectangle_points([-1.0, +1.0], [-1.0, +1.0])
    vals = quadpy.quadrilateral.integrate(eval_orthopolys, quad, scheme)
    # Put vals back into the tree structure:
    # len(approximate[k]) == k+1
    approximate = [
        vals[k * (k + 1) // 2 : (k + 1) * (k + 2) // 2]
        for k in range(scheme.degree + 2)
    ]

    exact = [numpy.zeros(k + 1) for k in range(scheme.degree + 2)]
    exact[0][0] = 2.0

    degree = check_degree_ortho(approximate, exact, abs_tol=tol)

    assert degree >= scheme.degree, "Observed: {}, expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [Product(quadpy.line_segment.GaussLegendre(5))])
def test_show(scheme):
    quadpy.quadrilateral.show(scheme)
    return


if __name__ == "__main__":
    # scheme_ = Product(quadpy.line_segment.GaussLegendre(6))
    scheme_ = quadpy.quadrilateral.HammerStroud("3-2")
    test_show(scheme_)
    test_scheme(scheme_, 1.0e-14)
