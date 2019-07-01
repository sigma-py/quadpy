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
    [
        quadpy.wedge.felippa_1(),
        quadpy.wedge.felippa_2(),
        quadpy.wedge.felippa_3(),
        quadpy.wedge.felippa_4(),
        quadpy.wedge.felippa_5(),
        quadpy.wedge.felippa_6(),
    ]
    + [
        quadpy.wedge.kubatko_yeager_maggi_1(),
        quadpy.wedge.kubatko_yeager_maggi_2a(),
        quadpy.wedge.kubatko_yeager_maggi_2b(),
        quadpy.wedge.kubatko_yeager_maggi_3a(),
        quadpy.wedge.kubatko_yeager_maggi_3b(),
        quadpy.wedge.kubatko_yeager_maggi_3c(),
        quadpy.wedge.kubatko_yeager_maggi_3d(),
        quadpy.wedge.kubatko_yeager_maggi_4a(),
        quadpy.wedge.kubatko_yeager_maggi_4b(),
        quadpy.wedge.kubatko_yeager_maggi_5a(),
        quadpy.wedge.kubatko_yeager_maggi_5b(),
        quadpy.wedge.kubatko_yeager_maggi_5c(),
        quadpy.wedge.kubatko_yeager_maggi_6a(),
        quadpy.wedge.kubatko_yeager_maggi_6b(),
        quadpy.wedge.kubatko_yeager_maggi_6c(),
        quadpy.wedge.kubatko_yeager_maggi_7a(),
        quadpy.wedge.kubatko_yeager_maggi_7b(),
        quadpy.wedge.kubatko_yeager_maggi_7c(),
        quadpy.wedge.kubatko_yeager_maggi_8a(),
        quadpy.wedge.kubatko_yeager_maggi_8b(),
        quadpy.wedge.kubatko_yeager_maggi_9(),
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
        lambda poly: scheme.integrate(poly, wedge),
        # lambda k: _integrate_exact(k, wedge),
        _integrate_monomial_over_unit_wedge,
        3,
        scheme.degree + 1,
    )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize("scheme", [quadpy.wedge.felippa_4()])
def test_show(scheme):
    scheme.show(backend="mpl")
    return


if __name__ == "__main__":
    scheme_ = quadpy.wedge.Felippa(2)
    test_scheme(scheme_)
    # test_show(scheme_)
    quadpy.wedge.show(scheme_, backend="vtk")
