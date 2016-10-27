# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3, check_degree
import numpy
import quadrature
import pytest
import sympy

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integrate_exact(k, wedge):
    def f(x):
        return x[0]**int(k[0]) * x[1]**int(k[1]) * x[2]**int(k[2])

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


@pytest.mark.parametrize(
    'scheme',
    [quadrature.wedge.Felippa(k) for k in range(1, 7)]
    )
def test_scheme(scheme):
    wedge = numpy.array([
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, +1.0],
        [1.0, 0.0, +1.0],
        [0.0, 1.0, +1.0],
        ])
    degree = check_degree(
            lambda poly: quadrature.wedge.integrate(
                poly, wedge, scheme
                ),
            lambda k: _integrate_exact(k, wedge),
            create_monomial_exponents3,
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.wedge.Felippa(4)]
    )
def test_show(scheme):
    wedge = numpy.array([
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, +1.0],
        [1.0, 0.0, +1.0],
        [0.0, 1.0, +1.0],
        ])
    quadrature.wedge.show(wedge, scheme)
    return


if __name__ == '__main__':
    scheme = quadrature.wedge.Felippa(2)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
