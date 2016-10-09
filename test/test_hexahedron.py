# -*- coding: utf-8 -*-
#
import numpy
import numpy.testing
import quadrature
import sympy

from test_triangle import _create_monomials

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
    det_J = sympy.simplify(sympy.det(J))
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.
    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))
    g_xi = sympy.simplify(f(pxi))
    exact = \
        sympy.integrate(
            sympy.integrate(
                sympy.integrate(abs_det_J * g_xi, (xi[2], -1, 1)),
                (xi[1], -1, 1)
            ),
            (xi[0], -1, 1)
        )
    return float(exact)


def test_generator():
    from quadrature.hexahedron import From1d
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
    schemes = [
        From1d(quadrature.line.Midpoint()),
        From1d(quadrature.line.Trapezoidal()),
        From1d(quadrature.line.GaussLegendre(1)),
        From1d(quadrature.line.GaussLegendre(2)),
        From1d(quadrature.line.NewtonCotesClosed(1)),
        From1d(quadrature.line.NewtonCotesClosed(2)),
        From1d(quadrature.line.NewtonCotesOpen(1)),
        From1d(quadrature.line.NewtonCotesOpen(2)),
        ]
    for scheme in schemes:
        yield check_scheme, scheme, hexa


def check_scheme(scheme, hexa):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for poly in _create_monomials(degree):
            exact_val = _integrate_exact(poly, hexa)
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


def test_show():
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
            # quadrature.line.Midpoint()
            # quadrature.line.Trapezoidal()
            # quadrature.line.NewtonCotesClosed(3)
            quadrature.line.NewtonCotesOpen(4)
            )
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
