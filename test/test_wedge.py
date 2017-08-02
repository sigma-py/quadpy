# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
import scipy.special

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
    '''Unit wedge given by the points
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, +1.0],
        [1.0, 0.0, +1.0],
        [0.0, 1.0, +1.0].
    '''
    if k[2] % 2 == 1:
        return 0.0
    return 2.0 * scipy.special.beta(k[0]+1, k[1]+2) / (k[1]+1) / (k[2]+1)


@pytest.mark.parametrize(
    'scheme',
    [quadpy.wedge.Felippa(k) for k in range(1, 7)]
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
            lambda poly: quadpy.wedge.integrate(
                poly, wedge, scheme
                ),
            # lambda k: _integrate_exact(k, wedge),
            _integrate_monomial_over_unit_wedge,
            lambda n: quadpy.helpers.partition(n, 3),
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.wedge.Felippa(4)]
    )
def test_show(scheme):
    quadpy.wedge.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.wedge.Felippa(6)
    test_scheme(scheme_)
    # test_show(scheme_)
