# -*- coding: utf-8 -*-
#
from __future__ import print_function

import pytest
import quadpy
from quadpy.quadrilateral import Product
import sympy

from helpers import check_degree


def _integrate_exact(f, quadrilateral):
    xi = sympy.DeferredVector('xi')
    pxi = quadrilateral[0] * 0.25*(1.0 + xi[0])*(1.0 + xi[1]) \
        + quadrilateral[1] * 0.25*(1.0 - xi[0])*(1.0 + xi[1]) \
        + quadrilateral[2] * 0.25*(1.0 - xi[0])*(1.0 - xi[1]) \
        + quadrilateral[3] * 0.25*(1.0 + xi[0])*(1.0 - xi[1])
    pxi = [
        sympy.expand(pxi[0]),
        sympy.expand(pxi[1]),
        ]
    # determinant of the transformation matrix
    det_J = \
        + sympy.diff(pxi[0], xi[0]) * sympy.diff(pxi[1], xi[1]) \
        - sympy.diff(pxi[1], xi[0]) * sympy.diff(pxi[0], xi[1])
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.
    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))

    g_xi = f(pxi)

    exact = sympy.integrate(
        sympy.integrate(abs_det_J * g_xi, (xi[1], -1, 1)),
        (xi[0], -1, 1)
        )
    return float(exact)


def _integrate_exact2(k, x0, x1, y0, y1):
    return 1.0/(k[0] + 1) * (x1**(k[0]+1) - x0**(k[0]+1)) \
        * 1.0/(k[1] + 1) * (y1**(k[1]+1) - y0**(k[1]+1))


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.quadrilateral.CoolsHaegemans1985(k), 1.0e-11)
        for k in range(1, 4)]
    + [(quadpy.quadrilateral.CoolsHaegemans1988(k), 1.0e-14) for k in [1, 2]]
    + [(quadpy.quadrilateral.Dunavant(k), 1.0e-14) for k in range(11)]
    + [(quadpy.quadrilateral.MorrowPatterson(k), 1.0e-7) for k in [1, 2]]
    + [(quadpy.quadrilateral.Stroud(k), 1.0e-14) for k in [
        'C2 1-1', 'C2 1-2',
        'C2 3-1', 'C2 3-2', 'C2 3-3', 'C2 3-4', 'C2 3-5',
        'C2 5-1', 'C2 5-2', 'C2 5-3', 'C2 5-4', 'C2 5-5', 'C2 5-6', 'C2 5-7'
        ]]
    + [(quadpy.quadrilateral.StroudN(k), 1.0e-14) for k in [
        'Cn 1-1', 'Cn 1-2',
        'Cn 2-1', 'Cn 2-2',
        'Cn 3-1', 'Cn 3-2', 'Cn 3-3', 'Cn 3-4', 'Cn 3-5', 'Cn 3-6',
        'Cn 5-2', 'Cn 5-3', 'Cn 5-4', 'Cn 5-5', 'Cn 5-6', 'Cn 5-7', 'Cn 5-9'
        ]]
    + [(quadpy.quadrilateral.StroudN(k), 1.0e-8) for k in ['Cn 7-1']]
    + [(quadpy.quadrilateral.WissmannBecker(k), 1.0e-14) for k in [
        '4-1', '4-2', '6-1', '6-2', '8-1', '8-2',
        ]]
    + [(Product(quadpy.line_segment.Midpoint()), 1.0e-14)]
    + [(Product(quadpy.line_segment.Trapezoidal()), 1.0e-14)]
    + [(Product(quadpy.line_segment.GaussLegendre(k)), 1.0e-14)
       for k in range(1, 5)
       ]
    + [(Product(quadpy.line_segment.NewtonCotesClosed(k)), 1.0e-14)
       for k in range(1, 5)
       ]
    + [(Product(quadpy.line_segment.NewtonCotesOpen(k)), 1.0e-14)
       for k in range(6)
       ]
    )
def test_scheme(scheme, tol, print_degree=False):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    x0 = -2.0
    x1 = +1.0
    y0 = -1.0
    y1 = +1.0
    quad = quadpy.quadrilateral.rectangle_points([-2.0, +1.0], [-1.0, +1.0])
    degree = check_degree(
            lambda poly: quadpy.quadrilateral.integrate(poly, quad, scheme),
            lambda k: _integrate_exact2(k, x0, x1, y0, y1),
            lambda n: quadpy.helpers.partition(n, 2),
            scheme.degree + 1,
            tol=tol
            )
    if print_degree:
        print('Detected degree {}, scheme degree {}.'.format(
            degree, scheme.degree
            ))
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [Product(quadpy.line_segment.GaussLegendre(5))]
    )
def test_show(scheme):
    quadpy.quadrilateral.show(scheme)
    return


if __name__ == '__main__':
    # scheme_ = Product(quadpy.line_segment.GaussLegendre(6))
    scheme_ = quadpy.quadrilateral.StroudN('Cn 7-1')
    print(scheme_.weights)
    print(scheme_.points)
    test_show(scheme_)
    test_scheme(scheme_, 1.0e-8, print_degree=True)
