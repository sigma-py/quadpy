import numpy
import pytest
import sympy

import orthopy
import quadpy
from helpers import check_degree_ortho

schemes = (
    [
        (quadpy.c2.albrecht_collatz_1(), 1.0e-14),
        (quadpy.c2.albrecht_collatz_2(), 1.0e-14),
        (quadpy.c2.albrecht_collatz_3(), 1.0e-14),
        (quadpy.c2.albrecht_collatz_4(), 1.0e-14),
    ]
    + [
        (quadpy.c2.cohen_gismalla_1(), 1.0e-6),
        (quadpy.c2.cohen_gismalla_2(), 1.0e-6),
    ]
    + [
        (quadpy.c2.cools_haegemans_1985_1(), 1.0e-10),
        (quadpy.c2.cools_haegemans_1985_2(), 1.0e-10),
        (quadpy.c2.cools_haegemans_1985_3(), 1.0e-10),
    ]
    + [
        (quadpy.c2.cools_haegemans_1988_1(), 1.0e-14),
        (quadpy.c2.cools_haegemans_1988_2(), 1.0e-14),
    ]
    + [
        (quadpy.c2.dunavant_00(), 1.0e-13),
        (quadpy.c2.dunavant_01(), 1.0e-13),
        (quadpy.c2.dunavant_02(), 1.0e-13),
        (quadpy.c2.dunavant_03(), 1.0e-13),
        (quadpy.c2.dunavant_04(), 1.0e-13),
        (quadpy.c2.dunavant_05(), 1.0e-13),
        (quadpy.c2.dunavant_06(), 1.0e-13),
        (quadpy.c2.dunavant_07(), 1.0e-13),
        (quadpy.c2.dunavant_08(), 1.0e-13),
        (quadpy.c2.dunavant_09(), 1.0e-13),
        (quadpy.c2.dunavant_10(), 1.0e-13),
    ]
    + [(quadpy.c2.franke_1(lmbda), 1.0e-13) for lmbda in [0.0, 1.0, -0.8]]
    + [
        (quadpy.c2.franke_2a(), 1.0e-13),
        (quadpy.c2.franke_2b(), 1.0e-13),
        (quadpy.c2.franke_3a(), 1.0e-13),
        (quadpy.c2.franke_3b(), 1.0e-13),
        (quadpy.c2.franke_3c(), 1.0e-13),
        (quadpy.c2.franke_5(), 1.0e-13),
        (quadpy.c2.franke_6(), 1.0e-13),
        (quadpy.c2.franke_8(), 1.0e-13),
    ]
    + [
        (quadpy.c2.hammer_stroud_1_2(), 1.0e-13),
        (quadpy.c2.hammer_stroud_2_2(), 1.0e-13),
        (quadpy.c2.hammer_stroud_3_2(), 1.0e-13),
    ]
    + [
        (quadpy.c2.morrow_patterson_1(), 1.0e-5),
        (quadpy.c2.morrow_patterson_2(), 1.0e-5),
    ]
    + [
        (quadpy.c2.stroud_c2_1_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_1_2(), 1.0e-13),
        (quadpy.c2.stroud_c2_3_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_3_2(), 1.0e-13),
        (quadpy.c2.stroud_c2_3_3(), 1.0e-13),
        (quadpy.c2.stroud_c2_3_4(), 1.0e-13),
        (quadpy.c2.stroud_c2_3_5(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_2(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_3(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_4(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_5(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_6(), 1.0e-13),
        (quadpy.c2.stroud_c2_5_7(), 1.0e-13),
        (quadpy.c2.stroud_c2_7_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_7_2(), 1.0e-13),
        (quadpy.c2.stroud_c2_7_3(), 1.0e-13),
        (quadpy.c2.stroud_c2_7_4(), 1.0e-13),
        (quadpy.c2.stroud_c2_7_5(), 1.0e-13),
        (quadpy.c2.stroud_c2_7_6(), 1.0e-13),
        (quadpy.c2.stroud_c2_9_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_11_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_11_2(), 1.0e-13),
        (quadpy.c2.stroud_c2_13_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_15_1(), 1.0e-13),
        (quadpy.c2.stroud_c2_15_2(), 1.0e-13),
    ]
    + [(quadpy.c2.haegemans_piessens(), 1.0e-14)]
    + [
        (quadpy.c2.piessens_haegemans_1(), 1.0e-14),
        (quadpy.c2.piessens_haegemans_2(), 1.0e-14),
    ]
    + [
        (quadpy.c2.schmid_2(), 1.0e-14),
        (quadpy.c2.schmid_4(), 1.0e-14),
        (quadpy.c2.schmid_6(), 1.0e-10),
    ]
    + [
        (quadpy.c2.sommariva_01(), 1.0e-10),
        (quadpy.c2.sommariva_02(), 1.0e-10),
        (quadpy.c2.sommariva_03(), 1.0e-10),
        (quadpy.c2.sommariva_04(), 1.0e-10),
        (quadpy.c2.sommariva_05(), 1.0e-10),
        (quadpy.c2.sommariva_06(), 1.0e-10),
        (quadpy.c2.sommariva_07(), 1.0e-10),
        (quadpy.c2.sommariva_08(), 1.0e-10),
        (quadpy.c2.sommariva_09(), 1.0e-10),
        (quadpy.c2.sommariva_10(), 1.0e-10),
        (quadpy.c2.sommariva_11(), 1.0e-10),
        (quadpy.c2.sommariva_12(), 1.0e-10),
        (quadpy.c2.sommariva_13(), 1.0e-10),
        (quadpy.c2.sommariva_14(), 1.0e-10),
        (quadpy.c2.sommariva_15(), 1.0e-10),
        (quadpy.c2.sommariva_16(), 1.0e-10),
        (quadpy.c2.sommariva_17(), 1.0e-10),
        (quadpy.c2.sommariva_18(), 1.0e-10),
        (quadpy.c2.sommariva_19(), 1.0e-10),
        (quadpy.c2.sommariva_20(), 1.0e-10),
        (quadpy.c2.sommariva_21(), 1.0e-10),
        (quadpy.c2.sommariva_22(), 1.0e-10),
        (quadpy.c2.sommariva_23(), 1.0e-10),
        (quadpy.c2.sommariva_24(), 1.0e-10),
        (quadpy.c2.sommariva_25(), 1.0e-10),
        (quadpy.c2.sommariva_26(), 1.0e-10),
        (quadpy.c2.sommariva_27(), 1.0e-10),
        (quadpy.c2.sommariva_28(), 1.0e-10),
        (quadpy.c2.sommariva_29(), 1.0e-10),
        (quadpy.c2.sommariva_30(), 1.0e-10),
        (quadpy.c2.sommariva_31(), 1.0e-10),
        (quadpy.c2.sommariva_32(), 1.0e-10),
        (quadpy.c2.sommariva_33(), 1.0e-10),
        (quadpy.c2.sommariva_34(), 1.0e-10),
        (quadpy.c2.sommariva_35(), 1.0e-10),
        (quadpy.c2.sommariva_36(), 1.0e-10),
        (quadpy.c2.sommariva_37(), 1.0e-10),
        (quadpy.c2.sommariva_38(), 1.0e-10),
        (quadpy.c2.sommariva_39(), 1.0e-10),
        (quadpy.c2.sommariva_40(), 1.0e-10),
        (quadpy.c2.sommariva_41(), 1.0e-10),
        (quadpy.c2.sommariva_42(), 1.0e-10),
        (quadpy.c2.sommariva_43(), 1.0e-10),
        (quadpy.c2.sommariva_44(), 1.0e-10),
        (quadpy.c2.sommariva_45(), 1.0e-10),
        (quadpy.c2.sommariva_46(), 1.0e-10),
        (quadpy.c2.sommariva_47(), 1.0e-10),
        (quadpy.c2.sommariva_48(), 1.0e-10),
        (quadpy.c2.sommariva_49(), 1.0e-10),
        (quadpy.c2.sommariva_50(), 1.0e-10),
        (quadpy.c2.sommariva_51(), 1.0e-10),
        (quadpy.c2.sommariva_52(), 1.0e-10),
        (quadpy.c2.sommariva_53(), 1.0e-10),
        (quadpy.c2.sommariva_54(), 1.0e-10),
        (quadpy.c2.sommariva_55(), 1.0e-10),
    ]
    + [(quadpy.c2.waldron(0.6, numpy.pi / 7), 1.0e-14)]
    + [
        (quadpy.c2.wissmann_becker_4_1(), 1.0e-14),
        (quadpy.c2.wissmann_becker_4_2(), 1.0e-14),
        (quadpy.c2.wissmann_becker_6_1(), 1.0e-14),
        (quadpy.c2.wissmann_becker_6_2(), 1.0e-14),
        (quadpy.c2.wissmann_becker_8_1(), 1.0e-14),
        (quadpy.c2.wissmann_becker_8_2(), 1.0e-14),
    ]
    + [
        (quadpy.c2.witherden_vincent_01(), 1.0e-14),
        (quadpy.c2.witherden_vincent_03(), 1.0e-14),
        (quadpy.c2.witherden_vincent_05(), 1.0e-14),
        # (quadpy.c2.witherden_vincent_07(), 1.0e-14), TODO
        (quadpy.c2.witherden_vincent_09(), 1.0e-14),
        (quadpy.c2.witherden_vincent_11(), 1.0e-14),
        (quadpy.c2.witherden_vincent_13(), 1.0e-14),
        (quadpy.c2.witherden_vincent_15(), 1.0e-14),
        (quadpy.c2.witherden_vincent_17(), 1.0e-14),
        (quadpy.c2.witherden_vincent_19(), 1.0e-14),
        (quadpy.c2.witherden_vincent_21(), 1.0e-14),
    ]
    + [(quadpy.c2.product(quadpy.c1.midpoint()), 1.0e-14)]
    + [(quadpy.c2.product(quadpy.c1.trapezoidal()), 1.0e-14)]
    + [
        (quadpy.c2.product(quadpy.c1.gauss_legendre(k)), 1.0e-14)
        for k in range(1, 5)
    ]
    + [
        (
            quadpy.c2.product(quadpy.c1.newton_cotes_closed(k)),
            1.0e-14,
        )
        for k in range(1, 5)
    ]
    + [
        (
            quadpy.c2.product(quadpy.c1.newton_cotes_open(k)),
            1.0e-14,
        )
        for k in range(6)
    ]
)


def _integrate_exact(f, c2):
    xi = sympy.DeferredVector("xi")
    pxi = (
        c2[0] * 0.25 * (1.0 + xi[0]) * (1.0 + xi[1])
        + c2[1] * 0.25 * (1.0 - xi[0]) * (1.0 + xi[1])
        + c2[2] * 0.25 * (1.0 - xi[0]) * (1.0 - xi[1])
        + c2[3] * 0.25 * (1.0 + xi[0]) * (1.0 - xi[1])
    )
    pxi = [sympy.expand(pxi[0]), sympy.expand(pxi[1])]
    # determinant of the transformation matrix
    det_J = +sympy.diff(pxi[0], xi[0]) * sympy.diff(pxi[1], xi[1]) - sympy.diff(
        pxi[1], xi[0]
    ) * sympy.diff(pxi[0], xi[1])
    # we cannot use abs(), see <https://github.com/sympy/sympy/issues/4212>.
    abs_det_J = sympy.Piecewise((det_J, det_J >= 0), (-det_J, det_J < 0))

    g_xi = f(pxi)

    exact = sympy.integrate(
        sympy.integrate(abs_det_J * g_xi, (xi[1], -1, 1)), (xi[0], -1, 1)
    )
    return float(exact)


def _integrate_exact2(k, x0, x1, y0, y1):
    return (
        1.0
        / (k[0] + 1)
        * (x1 ** (k[0] + 1) - x0 ** (k[0] + 1))
        * 1.0
        / (k[1] + 1)
        * (y1 ** (k[1] + 1) - y0 ** (k[1] + 1))
    )


@pytest.mark.parametrize("scheme,tol", schemes)
def test_scheme(scheme, tol):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    print(scheme.name)
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    def eval_orthopolys(x):
        return numpy.concatenate(
            orthopy.quadrilateral.tree(x, scheme.degree + 1, symbolic=False)
        )

    quad = quadpy.c2.rectangle_points([-1.0, +1.0], [-1.0, +1.0])
    vals = scheme.integrate(eval_orthopolys, quad)
    # Put vals back into the tree structure:
    # len(approximate[k]) == k+1
    approximate = [
        vals[k * (k + 1) // 2 : (k + 1) * (k + 2) // 2]
        for k in range(scheme.degree + 2)
    ]

    exact = [numpy.zeros(k + 1) for k in range(scheme.degree + 2)]
    exact[0][0] = 2.0

    degree = check_degree_ortho(approximate, exact, abs_tol=tol)

    assert degree >= scheme.degree, "Observed: {}, expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize(
    "scheme", [quadpy.c2.product(quadpy.c1.gauss_legendre(5))]
)
def test_show(scheme):
    scheme.show()
    return


if __name__ == "__main__":
    # scheme_ = Product(quadpy.c1.gauss_legendre(6))
    # scheme_ = quadpy.c2.HammerStroud("3-2")
    # scheme_ = quadpy.c2.Stroud["C2 3-2"]()
    # test_show(scheme_)
    # test_scheme(scheme_, 1.0e-14)
    from helpers import find_equal

    schemes_ = [scheme[0] for scheme in schemes]
    find_equal(schemes_)
