# -*- coding: utf-8 -*-
#
from __future__ import print_function

import numpy
import pytest
import sympy

import orthopy
import quadpy

from helpers import check_degree_ortho

schemes = (
    [
        (quadpy.quadrilateral.albrecht_collatz_1(), 1.0e-14),
        (quadpy.quadrilateral.albrecht_collatz_2(), 1.0e-14),
        (quadpy.quadrilateral.albrecht_collatz_3(), 1.0e-14),
        (quadpy.quadrilateral.albrecht_collatz_4(), 1.0e-14),
    ]
    + [
        (quadpy.quadrilateral.cohen_gismalla_1(), 1.0e-6),
        (quadpy.quadrilateral.cohen_gismalla_2(), 1.0e-6),
    ]
    + [
        (quadpy.quadrilateral.cools_haegemans_1985_1(), 1.0e-10),
        (quadpy.quadrilateral.cools_haegemans_1985_2(), 1.0e-10),
        (quadpy.quadrilateral.cools_haegemans_1985_3(), 1.0e-10),
    ]
    + [
        (quadpy.quadrilateral.cools_haegemans_1988_1(), 1.0e-14),
        (quadpy.quadrilateral.cools_haegemans_1988_2(), 1.0e-14),
    ]
    + [
        (quadpy.quadrilateral.dunavant_00(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_01(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_02(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_03(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_04(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_05(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_06(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_07(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_08(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_09(), 1.0e-13),
        (quadpy.quadrilateral.dunavant_10(), 1.0e-13),
    ]
    + [(quadpy.quadrilateral.franke_1(lmbda), 1.0e-13) for lmbda in [0.0, 1.0, -0.8]]
    + [
        (quadpy.quadrilateral.franke_2a(), 1.0e-13),
        (quadpy.quadrilateral.franke_3a(), 1.0e-13),
        (quadpy.quadrilateral.franke_3b(), 1.0e-13),
        (quadpy.quadrilateral.franke_3c(), 1.0e-13),
        (quadpy.quadrilateral.franke_5(), 1.0e-13),
        (quadpy.quadrilateral.franke_6(), 1.0e-13),
        (quadpy.quadrilateral.franke_8(), 1.0e-13),
    ]
    + [
        (quadpy.quadrilateral.hammer_stroud_1_2(), 1.0e-13),
        (quadpy.quadrilateral.hammer_stroud_2_2(), 1.0e-13),
        (quadpy.quadrilateral.hammer_stroud_3_2(), 1.0e-13),
    ]
    + [
        (quadpy.quadrilateral.morrow_patterson_1(), 1.0e-5),
        (quadpy.quadrilateral.morrow_patterson_2(), 1.0e-5),
    ]
    + [
        (quadpy.quadrilateral.stroud_c2_1_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_1_2(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_3_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_3_2(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_3_3(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_3_4(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_3_5(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_2(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_3(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_4(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_5(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_6(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_5_7(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_7_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_7_2(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_7_3(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_7_4(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_7_5(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_7_6(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_9_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_11_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_11_2(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_13_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_15_1(), 1.0e-13),
        (quadpy.quadrilateral.stroud_c2_15_2(), 1.0e-13),
    ]
    + [(quadpy.quadrilateral.haegemans_piessens(), 1.0e-14)]
    + [
        (quadpy.quadrilateral.piessens_haegemans_1(), 1.0e-14),
        (quadpy.quadrilateral.piessens_haegemans_2(), 1.0e-14),
    ]
    + [
        (quadpy.quadrilateral.schmid_2(), 1.0e-14),
        (quadpy.quadrilateral.schmid_4(), 1.0e-14),
        (quadpy.quadrilateral.schmid_6(), 1.0e-10),
    ]
    + [
        (quadpy.quadrilateral.sommariva_01(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_02(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_03(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_04(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_05(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_06(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_07(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_08(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_09(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_10(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_11(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_12(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_13(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_14(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_15(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_16(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_17(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_18(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_19(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_20(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_21(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_22(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_23(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_24(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_25(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_26(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_27(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_28(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_29(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_30(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_31(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_32(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_33(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_34(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_35(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_36(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_37(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_38(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_39(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_40(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_41(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_42(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_43(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_44(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_45(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_46(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_47(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_48(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_49(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_50(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_51(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_52(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_53(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_54(), 1.0e-10),
        (quadpy.quadrilateral.sommariva_55(), 1.0e-10),
    ]
    + [(quadpy.quadrilateral.waldron(0.6, numpy.pi / 7), 1.0e-14)]
    + [
        (quadpy.quadrilateral.wissmann_becker_4_1(), 1.0e-14),
        (quadpy.quadrilateral.wissmann_becker_4_2(), 1.0e-14),
        (quadpy.quadrilateral.wissmann_becker_6_1(), 1.0e-14),
        (quadpy.quadrilateral.wissmann_becker_6_2(), 1.0e-14),
        (quadpy.quadrilateral.wissmann_becker_8_1(), 1.0e-14),
        (quadpy.quadrilateral.wissmann_becker_8_2(), 1.0e-14),
    ]
    + [
        (quadpy.quadrilateral.witherden_vincent_01(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_03(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_05(), 1.0e-14),
        # (quadpy.quadrilateral.witherden_vincent_07(), 1.0e-14), TODO
        (quadpy.quadrilateral.witherden_vincent_09(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_11(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_13(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_15(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_17(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_19(), 1.0e-14),
        (quadpy.quadrilateral.witherden_vincent_21(), 1.0e-14),
    ]
    + [(quadpy.quadrilateral.product(quadpy.line_segment.Midpoint()), 1.0e-14)]
    + [(quadpy.quadrilateral.product(quadpy.line_segment.Trapezoidal()), 1.0e-14)]
    + [
        (quadpy.quadrilateral.product(quadpy.line_segment.GaussLegendre(k)), 1.0e-14)
        for k in range(1, 5)
    ]
    + [
        (
            quadpy.quadrilateral.product(quadpy.line_segment.NewtonCotesClosed(k)),
            1.0e-14,
        )
        for k in range(1, 5)
    ]
    + [
        (quadpy.quadrilateral.product(quadpy.line_segment.NewtonCotesOpen(k)), 1.0e-14)
        for k in range(6)
    ]
    + [
        # (quadpy.ncube.dobrodeev_1970(2), 1.0e-14),
        (quadpy.ncube.dobrodeev_1978(2), 1.0e-14),
        (quadpy.ncube.hammer_stroud_1n(2), 1.0e-14),
        (quadpy.ncube.hammer_stroud_2n(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_1_1(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_1_2(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_2_1(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_2_2(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_3_1(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_3_2(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_3_3(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_3_4(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_3_5(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_3_6(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_2(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_3(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_4(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_5(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_6(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_7(2), 1.0e-14),
        # (quadpy.ncube.stroud_cn_5_8(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_5_9(2), 1.0e-14),
        (quadpy.ncube.stroud_cn_7_1(2), 1.0e-14),
    ]
)


def _integrate_exact(f, quadrilateral):
    xi = sympy.DeferredVector("xi")
    pxi = (
        quadrilateral[0] * 0.25 * (1.0 + xi[0]) * (1.0 + xi[1])
        + quadrilateral[1] * 0.25 * (1.0 - xi[0]) * (1.0 + xi[1])
        + quadrilateral[2] * 0.25 * (1.0 - xi[0]) * (1.0 - xi[1])
        + quadrilateral[3] * 0.25 * (1.0 + xi[0]) * (1.0 - xi[1])
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

    quad = quadpy.quadrilateral.rectangle_points([-1.0, +1.0], [-1.0, +1.0])
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
    "scheme", [quadpy.quadrilateral.product(quadpy.line_segment.GaussLegendre(5))]
)
def test_show(scheme):
    scheme.show()
    return


if __name__ == "__main__":
    # scheme_ = Product(quadpy.line_segment.GaussLegendre(6))
    # scheme_ = quadpy.quadrilateral.HammerStroud("3-2")
    # scheme_ = quadpy.quadrilateral.Stroud["C2 3-2"]()
    # test_show(scheme_)
    # test_scheme(scheme_, 1.0e-14)
    from helpers import find_equal

    schemes_ = [scheme[0] for scheme in schemes]
    find_equal(schemes_)
