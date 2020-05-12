import matplotlib.pyplot as plt
import numpy
import pytest
import sympy

import orthopy
import quadpy
from helpers import check_degree_ortho

schemes = (
    [quadpy.c3.product(quadpy.c1.midpoint())]
    + [quadpy.c3.product(quadpy.c1.trapezoidal())]
    + [quadpy.c3.product(quadpy.c1.gauss_legendre(k)) for k in range(1, 6)]
    + [quadpy.c3.product(quadpy.c1.newton_cotes_closed(k)) for k in range(1, 5)]
    + [quadpy.c3.product(quadpy.c1.newton_cotes_open(k)) for k in range(5)]
    + [
        quadpy.c3.hammer_stroud_1_3(),
        quadpy.c3.hammer_stroud_2_3(),
        quadpy.c3.hammer_stroud_4_3(),
        quadpy.c3.hammer_stroud_5_3a(),
        quadpy.c3.hammer_stroud_5_3b(),
        quadpy.c3.hammer_stroud_6_3(),
        quadpy.c3.hammer_wymore(),
        quadpy.c3.mustard_lyness_blatt_1(),
        quadpy.c3.mustard_lyness_blatt_2(),
        quadpy.c3.mustard_lyness_blatt_3(),
        quadpy.c3.mustard_lyness_blatt_4(),
        quadpy.c3.mustard_lyness_blatt_5(),
        quadpy.c3.mustard_lyness_blatt_6(),
        quadpy.c3.mustard_lyness_blatt_7(),
        quadpy.c3.sadowsky(),
        quadpy.c3.stroud_c3_3_1(),
        quadpy.c3.stroud_c3_3_2(),
        quadpy.c3.stroud_c3_3_3(),
        quadpy.c3.stroud_c3_3_4(),
        quadpy.c3.stroud_c3_3_5(),
        quadpy.c3.stroud_c3_3_6(),
        quadpy.c3.stroud_c3_3_7(),
        quadpy.c3.stroud_c3_5_1(),
        quadpy.c3.stroud_c3_5_2(),
        quadpy.c3.stroud_c3_5_3(),
        quadpy.c3.stroud_c3_5_4(),
        quadpy.c3.stroud_c3_5_5(),
        quadpy.c3.stroud_c3_5_6(),
        quadpy.c3.stroud_c3_5_7(),
        quadpy.c3.stroud_c3_5_8(),
        quadpy.c3.stroud_c3_7_1a(),
        quadpy.c3.stroud_c3_7_1b(),
        quadpy.c3.stroud_c3_7_2(),
        quadpy.c3.stroud_c3_7_3(),
        quadpy.c3.stroud_1967(),
        quadpy.c3.tyler_1(),
        quadpy.c3.tyler_2(),
    ]
)


def _integrate_exact(f, hexa):
    xi = sympy.DeferredVector("xi")
    pxi = (
        +hexa[0] * 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
        + hexa[1] * 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
        + hexa[2] * 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
        + hexa[3] * 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
        + hexa[4] * 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
        + hexa[5] * 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
        + hexa[6] * 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])
        + hexa[7] * 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])
    )
    pxi = [sympy.expand(pxi[0]), sympy.expand(pxi[1]), sympy.expand(pxi[2])]
    # determinant of the transformation matrix
    J = sympy.Matrix(
        [
            [
                sympy.diff(pxi[0], xi[0]),
                sympy.diff(pxi[0], xi[1]),
                sympy.diff(pxi[0], xi[2]),
            ],
            [
                sympy.diff(pxi[1], xi[0]),
                sympy.diff(pxi[1], xi[1]),
                sympy.diff(pxi[1], xi[2]),
            ],
            [
                sympy.diff(pxi[2], xi[0]),
                sympy.diff(pxi[2], xi[1]),
                sympy.diff(pxi[2], xi[2]),
            ],
        ]
    )
    det_J = sympy.det(J)
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.

    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))
    g_xi = f(pxi)
    exact = sympy.integrate(
        sympy.integrate(
            sympy.integrate(abs_det_J * g_xi, (xi[2], -1, 1)), (xi[1], -1, 1)
        ),
        (xi[0], -1, 1),
    )
    return float(exact)


def _integrate_exact2(k, x0, x1, y0, y1, z0, z1):
    return (
        1.0
        / (k[0] + 1)
        * (x1 ** (k[0] + 1) - x0 ** (k[0] + 1))
        * 1.0
        / (k[1] + 1)
        * (y1 ** (k[1] + 1) - y0 ** (k[1] + 1))
        * 1.0
        / (k[2] + 1)
        * (z1 ** (k[2] + 1) - z0 ** (k[2] + 1))
    )


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme, tol=1.0e-14, print_degree=False):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    x = [-1.0, +1.0]
    y = [-1.0, +1.0]
    z = [-1.0, +1.0]
    hexa = quadpy.c3.cube_points(x, y, z)

    # degree = check_degree(
    #     lambda poly: scheme.integrate(poly, hexa),
    #     lambda k: _integrate_exact2(k, x[0], x[1], y[0], y[1], z[0], z[1]),
    #     3,
    #     scheme.degree + 1,
    #     tol=tol,
    # )
    # if print_degree:
    #     print("Detected degree {}, scheme degree {}.".format(degree, scheme.degree))
    # assert degree == scheme.degree, scheme.name

    def eval_orthopolys(x):
        return numpy.concatenate(
            orthopy.hexahedron.tree(x, scheme.degree + 1, symbolic=False)
        )

    vals = scheme.integrate(eval_orthopolys, hexa)
    # Put vals back into the tree structure:
    # len(approximate[k]) == k+1
    approximate = [
        vals[k * (k + 1) * (k + 2) // 6 : (k + 1) * (k + 2) * (k + 3) // 6]
        for k in range(scheme.degree + 2)
    ]

    exact = [numpy.zeros(len(s)) for s in approximate]
    exact[0][0] = numpy.sqrt(2.0) * 2

    degree = check_degree_ortho(approximate, exact, abs_tol=tol)

    assert degree >= scheme.degree, "{} -- Observed: {}, expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize(
    "scheme", [quadpy.c3.product(quadpy.c1.newton_cotes_closed(2))]
)
def test_show(scheme):
    scheme.show(backend="mpl")
    plt.close()
    return


if __name__ == "__main__":
    # scheme_ = Product(quadpy.c1.NewtonCotesOpen(5))
    # scheme_ = quadpy.c3.HammerStroud("6-3")
    # test_scheme(scheme_, 1.0e-14, print_degree=True)
    # test_show(scheme_)
    # scheme_.show(backend="vtk")
    from helpers import find_equal

    find_equal(schemes)
