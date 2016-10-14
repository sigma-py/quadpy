# -*- coding: utf-8 -*-
#
import numpy
import numpy.testing
import quadrature
import pytest
import sympy

from test_tetrahedron import _create_monomial_exponents

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integrate_exact(f, wedge):
    xi = sympy.DeferredVector('xi')
    pxi = \
        + wedge[0] * 0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2]) \
        + wedge[1] * 0.5 * xi[0] * (1.0-xi[2]) \
        + wedge[2] * 0.5 * xi[1] * (1.0-xi[2]) \
        + wedge[3] * 0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2]) \
        + wedge[4] * 0.5 * xi[0] * (1.0+xi[2]) \
        + wedge[5] * 0.5 * xi[1] * (1.0+xi[2])
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
                (xi[1], 0, 1-xi[0])
            ),
            (xi[0], 0, 1)
        )
    return float(exact)


@pytest.mark.parametrize('scheme', [
    quadrature.wedge.Felippa(1),
    quadrature.wedge.Felippa(2),
    quadrature.wedge.Felippa(3),
    quadrature.wedge.Felippa(4),
    quadrature.wedge.Felippa(5),
    quadrature.wedge.Felippa(6),
    ])
def test_scheme(scheme):
    wedge = numpy.array([
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, +1.0],
        [1.0, 0.0, +1.0],
        [0.0, 1.0, +1.0],
        ])
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in _create_monomial_exponents(degree):
            def poly(x):
                return x[0]**k[0] + x[1]**k[1] + x[2]**k[2]
            exact_val = _integrate_exact(poly, wedge)
            val = quadrature.wedge.integrate(
                    poly, wedge, scheme
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


# def test_show():
#     wedge = numpy.array([
#         [-1, -1, -1],
#         [+1, -1, -1],
#         [+1, +1, -1],
#         [-1, +1, -1],
#         [-1, -1, +1],
#         [+1, -1, +1],
#         [+1, +1, +1],
#         [-1, +1, +1],
#         ])
#     quadrature.hexahedron.show(
#         hexa,
#         quadrature.hexahedron.From1d(
#             # quadrature.line.Midpoint()
#             # quadrature.line.Trapezoidal()
#             # quadrature.line.NewtonCotesClosed(3)
#             quadrature.line.NewtonCotesOpen(4)
#             )
#         )
#     return


if __name__ == '__main__':
    # test_show()
    # plt.show()
    scheme = quadrature.wedge.Felippa(1)
    test_scheme(scheme)
