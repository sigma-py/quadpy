# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3
import numpy
import numpy.testing
import quadrature
import pytest
import sympy

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integrate_exact(f, pyra):
    # map the reference hexahedron [-1,1]^3 to the pyramid
    xi = sympy.DeferredVector('xi')
    pxi = \
        + pyra[0] * (1 - xi[0])*(1 - xi[1])*(1 - xi[2]) / 8 \
        + pyra[1] * (1 + xi[0])*(1 - xi[1])*(1 - xi[2]) / 8 \
        + pyra[2] * (1 + xi[0])*(1 + xi[1])*(1 - xi[2]) / 8 \
        + pyra[3] * (1 - xi[0])*(1 + xi[1])*(1 - xi[2]) / 8 \
        + pyra[4] * (1 + xi[2]) / 2

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
    # abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))
    # This is quite the leap of faith, but sympy will cowardly bail out
    # otherwise.
    abs_det_J = det_J

    exact = sympy.integrate(
        sympy.integrate(
            sympy.integrate(abs_det_J * f(pxi), (xi[2], -1, 1)),
            (xi[1], -1, +1)
            ),
        (xi[0], -1, +1)
        )

    return float(exact)


@pytest.mark.parametrize('scheme', [
    quadrature.pyramid.Felippa(1),
    quadrature.pyramid.Felippa(2),
    quadrature.pyramid.Felippa(3),
    quadrature.pyramid.Felippa(4),
    quadrature.pyramid.Felippa(5),
    quadrature.pyramid.Felippa(6),
    quadrature.pyramid.Felippa(7),
    quadrature.pyramid.Felippa(8),
    quadrature.pyramid.Felippa(9),
    ])
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    pyra = numpy.array([
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [0, 0, 1],
        ])
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in create_monomial_exponents3(degree):
            def poly(x):
                return x[0]**int(k[0]) * x[1]**int(k[1]) * x[2]**int(k[2])
            exact_val = _integrate_exact(poly, pyra)
            val = quadrature.pyramid.integrate(
                    poly, pyra, scheme
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
    pyra = numpy.array([
        [-1, -1, -1],
        [+1, -1, -1],
        [+1, +1, -1],
        [-1, +1, -1],
        [0, 0, 1],
        ])
    quadrature.pyramid.show(
        pyra,
        quadrature.pyramid.Felippa(9),
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
