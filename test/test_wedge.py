# -*- coding: utf-8 -*-
#
import numpy
import pytest
import scipy.special

import quadpy

from helpers import check_degree


# def _integrate_exact(k, wedge):
#     import sympy
#     def f(x):
#         return x[0]**int(k[0]) * x[1]**int(k[1]) * x[2]**int(k[2])
#
#     xi = sympy.DeferredVector('xi')
#     pxi = (
#         + wedge[0] * 0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2])
#         + wedge[1] * 0.5 * xi[0] * (1.0-xi[2])
#         + wedge[2] * 0.5 * xi[1] * (1.0-xi[2])
#         + wedge[3] * 0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2])
#         + wedge[4] * 0.5 * xi[0] * (1.0+xi[2])
#         + wedge[5] * 0.5 * xi[1] * (1.0+xi[2])
#         )
#     pxi = [
#         sympy.expand(pxi[0]),
#         sympy.expand(pxi[1]),
#         sympy.expand(pxi[2]),
#         ]
#     # determinant of the transformation matrix
#     J = sympy.Matrix([
#         [sympy.diff(pxi[0], xi[0]),
#          sympy.diff(pxi[0], xi[1]),
#          sympy.diff(pxi[0], xi[2])],
#         [sympy.diff(pxi[1], xi[0]),
#          sympy.diff(pxi[1], xi[1]),
#          sympy.diff(pxi[1], xi[2])],
#         [sympy.diff(pxi[2], xi[0]),
#          sympy.diff(pxi[2], xi[1]),
#          sympy.diff(pxi[2], xi[2])],
#         ])
#     det_J = sympy.det(J)
#     # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.
#     abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))
#     g_xi = f(pxi)
#     exact = \
#         sympy.integrate(
#             sympy.integrate(
#                 sympy.integrate(abs_det_J * g_xi, (xi[2], -1, 1)),
#                 (xi[1], 0, 1-xi[0])
#             ),
#             (xi[0], 0, 1)
#         )
#     return float(exact)


def _integrate_monomial_over_unit_wedge(k):
    """Unit wedge given by the points
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, +1.0],
        [1.0, 0.0, +1.0],
        [0.0, 1.0, +1.0].
    """
    if k[2] % 2 == 1:
        return 0.0
    return 2.0 * scipy.special.beta(k[0] + 1, k[1] + 2) / (k[1] + 1) / (k[2] + 1)


@pytest.mark.parametrize(
    "scheme",
    [scheme() for scheme in quadpy.wedge.Felippa.values()]
    + [
        quadpy.wedge.KubatkoYeagerMaggi(k)
        for k in [
            "1",
            "2a",
            "2b",
            "3a",
            "3b",
            "3c",
            "3d",
            "4a",
            "4b",
            "5a",
            "5b",
            "5c",
            "6a",
            "6b",
            "6c",
            "7a",
            "7b",
            "7c",
            "8a",
            "8b",
            "9",
        ]
    ],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    wedge = numpy.array(
        [
            [[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 1.0, -1.0]],
            [[0.0, 0.0, +1.0], [1.0, 0.0, +1.0], [0.0, 1.0, +1.0]],
        ]
    )

    degree = check_degree(
        lambda poly: quadpy.wedge.integrate(poly, wedge, scheme),
        # lambda k: _integrate_exact(k, wedge),
        _integrate_monomial_over_unit_wedge,
        3,
        scheme.degree + 1,
    )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize("scheme", [quadpy.wedge.Felippa[4]()])
def test_show(scheme):
    quadpy.wedge.show(scheme)
    return


if __name__ == "__main__":
    scheme_ = quadpy.wedge.Felippa(2)
    test_scheme(scheme_)
    # test_show(scheme_)
    quadpy.wedge.show(scheme_, backend="vtk")
