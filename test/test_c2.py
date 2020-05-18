import numpy
import pytest
import sympy

import orthopy
import quadpy
from helpers import check_degree_ortho

schemes = (
    [
        quadpy.c2.albrecht_collatz_1(),
        quadpy.c2.albrecht_collatz_2(),
        quadpy.c2.albrecht_collatz_3(),
        quadpy.c2.albrecht_collatz_4(),
    ]
    + [quadpy.c2.cohen_gismalla_1(), quadpy.c2.cohen_gismalla_2()]
    + [
        quadpy.c2.cools_haegemans_1985_1(),
        quadpy.c2.cools_haegemans_1985_2(),
        quadpy.c2.cools_haegemans_1985_3(),
    ]
    + [quadpy.c2.cools_haegemans_1988_1(), quadpy.c2.cools_haegemans_1988_2()]
    + [
        quadpy.c2.dunavant_00(),
        quadpy.c2.dunavant_01(),
        quadpy.c2.dunavant_02(),
        quadpy.c2.dunavant_03(),
        quadpy.c2.dunavant_04(),
        quadpy.c2.dunavant_05(),
        quadpy.c2.dunavant_06(),
        quadpy.c2.dunavant_07(),
        quadpy.c2.dunavant_08(),
        quadpy.c2.dunavant_09(),
        quadpy.c2.dunavant_10(),
    ]
    + [quadpy.c2.franke_1(lmbda) for lmbda in [0.0, 1.0, -0.8]]
    + [
        quadpy.c2.franke_2a(),
        quadpy.c2.franke_2b(),
        quadpy.c2.franke_3a(),
        quadpy.c2.franke_3b(),
        quadpy.c2.franke_3c(),
        quadpy.c2.franke_5(),
        quadpy.c2.franke_6(),
        quadpy.c2.franke_8(),
    ]
    + [
        quadpy.c2.hammer_stroud_1_2(),
        quadpy.c2.hammer_stroud_2_2(),
        quadpy.c2.hammer_stroud_3_2(),
    ]
    + [quadpy.c2.morrow_patterson_1(), quadpy.c2.morrow_patterson_2()]
    + [
        quadpy.c2.stroud_c2_1_1(),
        quadpy.c2.stroud_c2_1_2(),
        quadpy.c2.stroud_c2_3_1(),
        quadpy.c2.stroud_c2_3_2(),
        quadpy.c2.stroud_c2_3_3(),
        quadpy.c2.stroud_c2_3_4(),
        quadpy.c2.stroud_c2_3_5(),
        quadpy.c2.stroud_c2_5_1(),
        quadpy.c2.stroud_c2_5_2(),
        quadpy.c2.stroud_c2_5_3(),
        quadpy.c2.stroud_c2_5_4(),
        quadpy.c2.stroud_c2_5_5(),
        quadpy.c2.stroud_c2_5_6(),
        quadpy.c2.stroud_c2_5_7(),
        quadpy.c2.stroud_c2_7_1(),
        quadpy.c2.stroud_c2_7_2(),
        quadpy.c2.stroud_c2_7_3(),
        quadpy.c2.stroud_c2_7_4(),
        quadpy.c2.stroud_c2_7_5(),
        quadpy.c2.stroud_c2_7_6(),
        quadpy.c2.stroud_c2_9_1(),
        quadpy.c2.stroud_c2_11_1(),
        quadpy.c2.stroud_c2_11_2(),
        quadpy.c2.stroud_c2_13_1(),
        quadpy.c2.stroud_c2_15_1(),
        quadpy.c2.stroud_c2_15_2(),
    ]
    + [quadpy.c2.haegemans_piessens()]
    + [quadpy.c2.piessens_haegemans_1(), quadpy.c2.piessens_haegemans_2()]
    + [quadpy.c2.schmid_2(), quadpy.c2.schmid_4(), quadpy.c2.schmid_6()]
    + [
        quadpy.c2.sommariva_01(),
        quadpy.c2.sommariva_02(),
        quadpy.c2.sommariva_03(),
        quadpy.c2.sommariva_04(),
        quadpy.c2.sommariva_05(),
        quadpy.c2.sommariva_06(),
        quadpy.c2.sommariva_07(),
        quadpy.c2.sommariva_08(),
        quadpy.c2.sommariva_09(),
        quadpy.c2.sommariva_10(),
        quadpy.c2.sommariva_11(),
        quadpy.c2.sommariva_12(),
        quadpy.c2.sommariva_13(),
        quadpy.c2.sommariva_14(),
        quadpy.c2.sommariva_15(),
        quadpy.c2.sommariva_16(),
        quadpy.c2.sommariva_17(),
        quadpy.c2.sommariva_18(),
        quadpy.c2.sommariva_19(),
        quadpy.c2.sommariva_20(),
        quadpy.c2.sommariva_21(),
        quadpy.c2.sommariva_22(),
        quadpy.c2.sommariva_23(),
        quadpy.c2.sommariva_24(),
        quadpy.c2.sommariva_25(),
        quadpy.c2.sommariva_26(),
        quadpy.c2.sommariva_27(),
        quadpy.c2.sommariva_28(),
        quadpy.c2.sommariva_29(),
        quadpy.c2.sommariva_30(),
        quadpy.c2.sommariva_31(),
        quadpy.c2.sommariva_32(),
        quadpy.c2.sommariva_33(),
        quadpy.c2.sommariva_34(),
        quadpy.c2.sommariva_35(),
        quadpy.c2.sommariva_36(),
        quadpy.c2.sommariva_37(),
        quadpy.c2.sommariva_38(),
        quadpy.c2.sommariva_39(),
        quadpy.c2.sommariva_40(),
        quadpy.c2.sommariva_41(),
        quadpy.c2.sommariva_42(),
        quadpy.c2.sommariva_43(),
        quadpy.c2.sommariva_44(),
        quadpy.c2.sommariva_45(),
        quadpy.c2.sommariva_46(),
        quadpy.c2.sommariva_47(),
        quadpy.c2.sommariva_48(),
        quadpy.c2.sommariva_49(),
        quadpy.c2.sommariva_50(),
        quadpy.c2.sommariva_51(),
        quadpy.c2.sommariva_52(),
        quadpy.c2.sommariva_53(),
        quadpy.c2.sommariva_54(),
        quadpy.c2.sommariva_55(),
    ]
    + [quadpy.c2.waldron(0.6, numpy.pi / 7)]
    + [
        quadpy.c2.wissmann_becker_4_1(),
        quadpy.c2.wissmann_becker_4_2(),
        quadpy.c2.wissmann_becker_6_1(),
        quadpy.c2.wissmann_becker_6_2(),
        quadpy.c2.wissmann_becker_8_1(),
        quadpy.c2.wissmann_becker_8_2(),
    ]
    + [
        quadpy.c2.witherden_vincent_01(),
        quadpy.c2.witherden_vincent_03(),
        quadpy.c2.witherden_vincent_05(),
        # (quadpy.c2.witherden_vincent_07(), 1.0e-14), TODO
        quadpy.c2.witherden_vincent_09(),
        quadpy.c2.witherden_vincent_11(),
        quadpy.c2.witherden_vincent_13(),
        quadpy.c2.witherden_vincent_15(),
        quadpy.c2.witherden_vincent_17(),
        quadpy.c2.witherden_vincent_19(),
        quadpy.c2.witherden_vincent_21(),
    ]
    + [quadpy.c2.product(quadpy.c1.midpoint())]
    + [quadpy.c2.product(quadpy.c1.trapezoidal())]
    + [quadpy.c2.product(quadpy.c1.gauss_legendre(k)) for k in range(1, 5)]
    + [quadpy.c2.product(quadpy.c1.newton_cotes_closed(k)) for k in range(1, 5)]
    + [quadpy.c2.product(quadpy.c1.newton_cotes_open(k)) for k in range(6)]
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


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no longer be
    # integrated exactly. The scheme's degree is `d-1`.
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

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

    degree, err = check_degree_ortho(approximate, exact, abs_tol=scheme.test_tolerance)

    assert (
        degree >= scheme.degree
    ), "{} -- Observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.c2.product(quadpy.c1.gauss_legendre(5))])
def test_show(scheme):
    scheme.show()


if __name__ == "__main__":
    # scheme_ = Product(quadpy.c1.gauss_legendre(6))
    # scheme_ = quadpy.c2.HammerStroud("3-2")
    # scheme_ = quadpy.c2.Stroud["C2 3-2"]()
    # test_show(scheme_)
    # test_scheme(scheme_, 1.0e-14)
    from helpers import find_equal

    schemes_ = [scheme[0] for scheme in schemes]
    find_equal(schemes_)
