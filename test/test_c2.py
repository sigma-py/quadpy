import numpy
import orthopy
import pytest

import quadpy

schemes = (
    [
        quadpy.c2.albrecht_collatz_1(),
        quadpy.c2.albrecht_collatz_2(),
        quadpy.c2.albrecht_collatz_3(),
        quadpy.c2.albrecht_collatz_4(),
    ]
    + [quadpy.c2.cohen_gismalla_1(), quadpy.c2.cohen_gismalla_2()]
    + [
        quadpy.c2.cools_haegemans_1985_9_1(),
        quadpy.c2.cools_haegemans_1985_13_1(),
        quadpy.c2.cools_haegemans_1985_13_2(),
        quadpy.c2.cools_haegemans_1985_13_3(),
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
        quadpy.c2.franke_7(),
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
        quadpy.c2.witherden_vincent_07(),
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


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    quad = quadpy.c2.rectangle_points([-1.0, +1.0], [-1.0, +1.0])

    evaluator = orthopy.cn.Eval(scheme.points.T)

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), quad)
        exact = evaluator.int_p0 * 4 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        max_err = max(max_err, numpy.max(err))
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator)[0], quad)
            exact = 2.0 if i == 0 else 0.0
            err = numpy.abs(approximate - exact)
            max_err = max(max_err, numpy.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
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
