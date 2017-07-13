# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents2, check_degree

from matplotlib import pyplot as plt
import numpy
import pytest
import quadpy
from quadpy.quadrilateral import Product
import sympy


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
    'scheme',
    [Product(quadpy.line_segment.Midpoint())]
    + [Product(quadpy.line_segment.Trapezoidal())]
    + [quadpy.quadrilateral.Stroud(k) for k in range(1, 7)]
    + [Product(quadpy.line_segment.GaussLegendre(k)) for k in range(1, 5)]
    + [Product(quadpy.line_segment.NewtonCotesClosed(k))
        for k in range(1, 5)
       ]
    + [Product(quadpy.line_segment.NewtonCotesOpen(k)) for k in range(6)]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    x0 = -2.0
    x1 = +1.0
    y0 = -1.0
    y1 = +1.0
    quadrilateral = numpy.array([
        [x0, y0],
        [x1, y0],
        [x1, y1],
        [x0, y1],
        ])
    degree = check_degree(
            lambda poly: quadpy.quadrilateral.integrate(
                poly, quadrilateral, scheme
                ),
            lambda k: _integrate_exact2(k, x0, x1, y0, y1),
            create_monomial_exponents2,
            scheme.degree + 1
            )
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
    scheme_ = Product(quadpy.line_segment.GaussLegendre(6))
    # scheme_ = quadpy.quadrilateral.Stroud(6)
    test_scheme(scheme_)
    test_show(scheme_)
    plt.show()
