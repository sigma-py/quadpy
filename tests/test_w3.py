import numpy as np
import pytest
import scipy.special
from helpers import check_degree

import quadpy

# def _integrate_exact(k, w3):
#     import sympy
#     def f(x):
#         return x[0]**int(k[0]) * x[1]**int(k[1]) * x[2]**int(k[2])
#
#     xi = sympy.DeferredVector('xi')
#     pxi = (
#         + w3[0] * 0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2])
#         + w3[1] * 0.5 * xi[0] * (1.0-xi[2])
#         + w3[2] * 0.5 * xi[1] * (1.0-xi[2])
#         + w3[3] * 0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2])
#         + w3[4] * 0.5 * xi[0] * (1.0+xi[2])
#         + w3[5] * 0.5 * xi[1] * (1.0+xi[2])
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


def _integrate_monomial_over_unit_w3(k):
    """Unit w3 given by the points
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
        quadpy.w3.felippa_1(),
        quadpy.w3.felippa_2(),
        quadpy.w3.felippa_3(),
        quadpy.w3.felippa_4(),
        quadpy.w3.felippa_5(),
        quadpy.w3.felippa_6(),
    ]
    + [
        quadpy.w3.kubatko_yeager_maggi_1(),
        quadpy.w3.kubatko_yeager_maggi_2a(),
        quadpy.w3.kubatko_yeager_maggi_2b(),
        quadpy.w3.kubatko_yeager_maggi_3a(),
        quadpy.w3.kubatko_yeager_maggi_3b(),
        quadpy.w3.kubatko_yeager_maggi_3c(),
        quadpy.w3.kubatko_yeager_maggi_3d(),
        quadpy.w3.kubatko_yeager_maggi_4a(),
        quadpy.w3.kubatko_yeager_maggi_4b(),
        quadpy.w3.kubatko_yeager_maggi_5a(),
        quadpy.w3.kubatko_yeager_maggi_5b(),
        quadpy.w3.kubatko_yeager_maggi_5c(),
        quadpy.w3.kubatko_yeager_maggi_6a(),
        quadpy.w3.kubatko_yeager_maggi_6b(),
        quadpy.w3.kubatko_yeager_maggi_6c(),
        quadpy.w3.kubatko_yeager_maggi_7a(),
        quadpy.w3.kubatko_yeager_maggi_7b(),
        quadpy.w3.kubatko_yeager_maggi_7c(),
        quadpy.w3.kubatko_yeager_maggi_8a(),
        quadpy.w3.kubatko_yeager_maggi_8b(),
        quadpy.w3.kubatko_yeager_maggi_9(),
    ],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    print(scheme)

    w3 = np.array(
        [
            [[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 1.0, -1.0]],
            [[0.0, 0.0, +1.0], [1.0, 0.0, +1.0], [0.0, 1.0, +1.0]],
        ]
    )

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, w3),
        # lambda k: _integrate_exact(k, w3),
        _integrate_monomial_over_unit_w3,
        3,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{}  --  observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.w3.felippa_4()])
def test_show(scheme):
    plt = scheme.show(backend="mpl")
    plt.close()


if __name__ == "__main__":
    scheme_ = quadpy.w3.Felippa(2)
    test_scheme(scheme_)
    # test_show(scheme_)
    quadpy.w3.show(scheme_, backend="vtk")
