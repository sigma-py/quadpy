# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3
import numpy
import numpy.testing
import quadrature
from quadrature.hexahedron import From1d
import pytest
import sympy

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


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


def _integrate_exact2(k, x0, x1, y0, y1, z0, z1):
    return 1.0/(k[0] + 1) * (x1**(k[0]+1) - x0**(k[0]+1)) \
        * 1.0/(k[1] + 1) * (y1**(k[1]+1) - y0**(k[1]+1)) \
        * 1.0/(k[2] + 1) * (z1**(k[2]+1) - z0**(k[2]+1))


@pytest.mark.parametrize(
    'scheme',
    [From1d(quadrature.line_segment.Midpoint())]
    + [From1d(quadrature.line_segment.Trapezoidal())]
    + [From1d(quadrature.line_segment.GaussLegendre(k)) for k in range(1, 6)]
    + [From1d(quadrature.line_segment.NewtonCotesClosed(k))
        for k in range(1, 5)
       ]
    + [From1d(quadrature.line_segment.NewtonCotesOpen(k)) for k in range(5)]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
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
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in create_monomial_exponents3(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1] * x[2]**k[2]
            # exact_val = _integrate_exact(poly, hexa)
            exact_val = _integrate_exact2(k, x0, x1, y0, y1, z0, z1)
            val = quadrature.hexahedron.integrate(
                    poly, hexa, scheme
                    )
            if abs(exact_val - val) > 1.0e-10:
                success = False
                break
        if not success:
            break
        if degree >= max_degree:
            break
        degree += 1
    numpy.testing.assert_equal(degree-1, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [From1d(quadrature.line_segment.NewtonCotesClosed(3))]
    )
def test_show(scheme):
    hexa = numpy.array([
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [-1, -1, +1],
        [+1, -1, +1],
        [+1, +1, +1],
        [-1, +1, +1],
        ])
    quadrature.hexahedron.show(
        hexa,
        quadrature.hexahedron.From1d(
            # quadrature.line_segment.Midpoint()
            # quadrature.line_segment.Trapezoidal()
            quadrature.line_segment.NewtonCotesClosed(2)
            # quadrature.line_segment.NewtonCotesOpen(2)
            )
        )
    return


if __name__ == '__main__':
    # scheme = From1d(quadrature.line_segment.NewtonCotesOpen(2))
    scheme = From1d(quadrature.line_segment.GaussLegendre(5))
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
