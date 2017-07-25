# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3, check_degree
from matplotlib import pyplot as plt
import numpy
import quadpy
from quadpy.hexahedron import Product
import pytest
import sympy


def _integrate_exact(f, hexa):
    xi = sympy.DeferredVector('xi')
    pxi = \
        + hexa[0] * 0.125*(1.0 - xi[0])*(1.0 - xi[1])*(1.0 - xi[2]) \
        + hexa[1] * 0.125*(1.0 + xi[0])*(1.0 - xi[1])*(1.0 - xi[2]) \
        + hexa[2] * 0.125*(1.0 + xi[0])*(1.0 + xi[1])*(1.0 - xi[2]) \
        + hexa[3] * 0.125*(1.0 - xi[0])*(1.0 + xi[1])*(1.0 - xi[2]) \
        + hexa[4] * 0.125*(1.0 - xi[0])*(1.0 - xi[1])*(1.0 + xi[2]) \
        + hexa[5] * 0.125*(1.0 + xi[0])*(1.0 - xi[1])*(1.0 + xi[2]) \
        + hexa[6] * 0.125*(1.0 + xi[0])*(1.0 + xi[1])*(1.0 + xi[2]) \
        + hexa[7] * 0.125*(1.0 - xi[0])*(1.0 + xi[1])*(1.0 + xi[2])
    pxi = [
        sympy.expand(pxi[0]),
        sympy.expand(pxi[1]),
        sympy.expand(pxi[2]),
        ]
    # determinant of the transformation matrix
    J = sympy.Matrix([
        [sympy.diff(pxi[0], xi[0]),
         sympy.diff(pxi[0], xi[1]),
         sympy.diff(pxi[0], xi[2])],
        [sympy.diff(pxi[1], xi[0]),
         sympy.diff(pxi[1], xi[1]),
         sympy.diff(pxi[1], xi[2])],
        [sympy.diff(pxi[2], xi[0]),
         sympy.diff(pxi[2], xi[1]),
         sympy.diff(pxi[2], xi[2])],
        ])
    det_J = sympy.det(J)
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.
    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))
    g_xi = f(pxi)
    exact = \
        sympy.integrate(
            sympy.integrate(
                sympy.integrate(abs_det_J * g_xi, (xi[2], -1, 1)),
                (xi[1], -1, 1)
            ),
            (xi[0], -1, 1)
        )
    return float(exact)


# pylint: disable=too-many-arguments
def _integrate_exact2(k, x0, x1, y0, y1, z0, z1):
    return 1.0/(k[0] + 1) * (x1**(k[0]+1) - x0**(k[0]+1)) \
        * 1.0/(k[1] + 1) * (y1**(k[1]+1) - y0**(k[1]+1)) \
        * 1.0/(k[2] + 1) * (z1**(k[2]+1) - z0**(k[2]+1))


@pytest.mark.parametrize(
    'scheme',
    [Product(quadpy.line_segment.Midpoint())]
    + [Product(quadpy.line_segment.Trapezoidal())]
    + [Product(quadpy.line_segment.GaussLegendre(k)) for k in range(1, 6)]
    + [Product(quadpy.line_segment.NewtonCotesClosed(k))
        for k in range(1, 5)
       ]
    + [Product(quadpy.line_segment.NewtonCotesOpen(k)) for k in range(5)]
    )
def test_scheme(scheme):
    x0 = -1
    x1 = +1
    y0 = -1
    y1 = +1
    z0 = -1
    z1 = +1
    hexa = numpy.array([
        [x0, y0, z0],
        [x1, y0, z0],
        [x1, y1, z0],
        [x0, y1, z0],
        [x0, y0, z1],
        [x1, y0, z1],
        [x1, y1, z1],
        [x0, y1, z1],
        ])
    degree = check_degree(
            lambda poly: quadpy.hexahedron.integrate(poly, hexa, scheme),
            lambda k: _integrate_exact2(k, x0, x1, y0, y1, z0, z1),
            create_monomial_exponents3,
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [Product(quadpy.line_segment.NewtonCotesClosed(3))]
    )
def test_show(scheme):
    quadpy.hexahedron.show(scheme)
    return


if __name__ == '__main__':
    # scheme_ = Product(quadpy.line_segment.NewtonCotesOpen(2))
    scheme_ = quadpy.hexahedron.StroudN('Cn 5-3')
    print(scheme_.weights)
    print(scheme_.points)
    test_scheme(scheme_)
    test_show(scheme_)
    plt.show()
