# -*- coding: utf-8 -*-
#
import numpy
import pytest

# from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
import sympy

import orthopy
import quadpy

from helpers import check_degree_ortho

schemes_tol = [
    (quadpy.triangle.albrecht_collatz(), 1.0e-14),
    (quadpy.triangle.centroid(), 1.0e-14),
    (quadpy.triangle.cools_haegemans_1(), 1.0e-13),
    # (quadpy.triangle.cools_haegemans_2(), 1.0e-13),
    (quadpy.triangle.cubtri(), 1.0e-14),
    (quadpy.triangle.berntsen_espelid_1(), 1.0e-11),
    (quadpy.triangle.berntsen_espelid_2(), 1.0e-11),
    (quadpy.triangle.berntsen_espelid_3(), 1.0e-11),
    (quadpy.triangle.berntsen_espelid_4(), 1.0e-11),
    (quadpy.triangle.dcutri(), 1.0e-11),
    (quadpy.triangle.dunavant_01(), 1.0e-14),
    (quadpy.triangle.dunavant_02(), 1.0e-14),
    (quadpy.triangle.dunavant_03(), 1.0e-14),
    (quadpy.triangle.dunavant_04(), 1.0e-14),
    (quadpy.triangle.dunavant_05(), 1.0e-12),
    (quadpy.triangle.dunavant_06(), 1.0e-14),
    (quadpy.triangle.dunavant_07(), 1.0e-14),
    (quadpy.triangle.dunavant_08(), 1.0e-14),
    (quadpy.triangle.dunavant_09(), 1.0e-14),
    (quadpy.triangle.dunavant_10(), 1.0e-13),
    (quadpy.triangle.dunavant_11(), 1.0e-12),
    (quadpy.triangle.dunavant_12(), 1.0e-14),
    (quadpy.triangle.dunavant_13(), 1.0e-14),
    (quadpy.triangle.dunavant_14(), 1.0e-14),
    (quadpy.triangle.dunavant_15(), 1.0e-13),
    (quadpy.triangle.dunavant_16(), 1.0e-13),
    (quadpy.triangle.dunavant_17(), 1.0e-13),
    (quadpy.triangle.dunavant_18(), 1.0e-12),
    (quadpy.triangle.dunavant_19(), 1.0e-13),
    (quadpy.triangle.dunavant_20(), 1.0e-13),
    (quadpy.triangle.franke_09(), 1.0e-14),
    (quadpy.triangle.franke_10(), 1.0e-14),
    (quadpy.triangle.gatermann(), 1.0e-12),
    (quadpy.triangle.griener_schmid_1(), 1.0e-14),
    (quadpy.triangle.griener_schmid_2(), 1.0e-14),
    (quadpy.triangle.hammer_marlowe_stroud_1(), 1.0e-14),
    (quadpy.triangle.hammer_marlowe_stroud_2(), 1.0e-14),
    (quadpy.triangle.hammer_marlowe_stroud_3(), 1.0e-14),
    (quadpy.triangle.hammer_marlowe_stroud_4(), 1.0e-14),
    (quadpy.triangle.hammer_marlowe_stroud_5(), 1.0e-14),
    (quadpy.triangle.hillion_01(), 1.0e-14),
    (quadpy.triangle.hillion_02(), 1.0e-14),
    (quadpy.triangle.hillion_03(), 1.0e-14),
    (quadpy.triangle.hillion_04(), 1.0e-14),
    (quadpy.triangle.hillion_05(), 1.0e-14),
    (quadpy.triangle.hillion_06(), 1.0e-14),
    (quadpy.triangle.hillion_07(), 1.0e-14),
    (quadpy.triangle.hillion_08(), 1.0e-14),
    (quadpy.triangle.hillion_09(), 1.0e-14),
    (quadpy.triangle.hillion_10(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_01(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_02a(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_02b(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_03(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_04(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_05(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_06(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_07(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_08(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_09(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_10(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_11(), 1.0e-13),
    (quadpy.triangle.laursen_gellert_12(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_13(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_14(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_15a(), 1.0e-14),
    (quadpy.triangle.laursen_gellert_15b(), 1.0e-14),
    (quadpy.triangle.lether(1), 1.0e-14),
    (quadpy.triangle.lether(2), 1.0e-14),
    (quadpy.triangle.lether(3), 1.0e-14),
    (quadpy.triangle.lether(4), 1.0e-14),
    (quadpy.triangle.lether(5), 1.0e-14),
    (quadpy.triangle.lether(6), 1.0e-14),
    (quadpy.triangle.lether(7), 1.0e-14),
    (quadpy.triangle.lether(8), 1.0e-14),
    (quadpy.triangle.lether(9), 1.0e-14),
    (quadpy.triangle.liu_vinokur_01(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_02(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_03(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_04(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_05(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_06(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_07(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_08(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_09(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_10(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_11(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_12(), 1.0e-14),
    (quadpy.triangle.liu_vinokur_13(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_01(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_02(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_03(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_04(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_05(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_06(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_07(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_08(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_09(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_10(), 1.0e-12),
    (quadpy.triangle.lyness_jespersen_11(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_12(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_13(), 1.0e-13),
    (quadpy.triangle.lyness_jespersen_14(), 1.0e-12),
    (quadpy.triangle.lyness_jespersen_15(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_16(), 1.0e-11),
    (quadpy.triangle.lyness_jespersen_17(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_18(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_19(), 1.0e-14),
    (quadpy.triangle.lyness_jespersen_20(), 1.0e-12),
    (quadpy.triangle.lyness_jespersen_21(), 1.0e-13),
    (quadpy.triangle.newton_cotes_closed(1), 1.0e-14),
    (quadpy.triangle.newton_cotes_closed(2), 1.0e-14),
    (quadpy.triangle.newton_cotes_closed(3), 1.0e-14),
    (quadpy.triangle.newton_cotes_closed(4), 1.0e-14),
    (quadpy.triangle.newton_cotes_closed(5), 1.0e-14),
    (quadpy.triangle.newton_cotes_closed(6), 1.0e-13),
    (quadpy.triangle.newton_cotes_open(1), 1.0e-14),
    (quadpy.triangle.newton_cotes_open(2), 1.0e-14),
    (quadpy.triangle.newton_cotes_open(3), 1.0e-14),
    (quadpy.triangle.newton_cotes_open(4), 1.0e-14),
    (quadpy.triangle.newton_cotes_open(5), 1.0e-13),
    (quadpy.triangle.newton_cotes_open(6), 1.0e-12),
    (quadpy.triangle.papanicolopulos_sym_0(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_1(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_2(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_3(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_4(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_5(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_6(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_7(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_sym_8(), 1.0e-13),
    (quadpy.triangle.papanicolopulos_rot_08(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_09(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_10(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_11(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_12(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_13(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_14(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_15(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_16(), 1.0e-14),
    (quadpy.triangle.papanicolopulos_rot_17(), 1.0e-14),
    (quadpy.triangle.seven_point(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_01(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_02(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_03(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_04(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_05(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_06(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_07(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_08(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_09(), 1.0e-14),
    (quadpy.triangle.strang_fix_cowper_10(), 1.0e-14),
    (quadpy.triangle.stroud_t2_3_1(), 1.0e-14),
    (quadpy.triangle.stroud_t2_5_1(), 1.0e-14),
    (quadpy.triangle.stroud_t2_7_1(), 1.0e-12),
    (quadpy.triangle.taylor_wingate_bos_1(), 1.0e-14),
    (quadpy.triangle.taylor_wingate_bos_2(), 1.0e-12),
    (quadpy.triangle.taylor_wingate_bos_4(), 1.0e-12),
    (quadpy.triangle.taylor_wingate_bos_5(), 1.0e-12),
    (quadpy.triangle.taylor_wingate_bos_8(), 1.0e-12),
    (quadpy.triangle.vertex(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_00(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_01(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_02(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_03(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_04(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_05(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_06(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_07(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_08(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_09(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_10(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_11(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_12(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_13(), 1.0e-13),
    (quadpy.triangle.vioreanu_rokhlin_14(), 1.0e-13),
    (quadpy.triangle.vioreanu_rokhlin_15(), 1.0e-13),
    (quadpy.triangle.vioreanu_rokhlin_16(), 1.0e-12),
    (quadpy.triangle.vioreanu_rokhlin_17(), 1.0e-14),
    (quadpy.triangle.vioreanu_rokhlin_18(), 1.0e-13),
    (quadpy.triangle.vioreanu_rokhlin_19(), 1.0e-11),
    (quadpy.triangle.walkington_p5(), 1.0e-14),
    (quadpy.triangle.wandzura_xiao_1(), 1.0e-14),
    (quadpy.triangle.wandzura_xiao_2(), 1.0e-14),
    (quadpy.triangle.wandzura_xiao_3(), 1.0e-14),
    (quadpy.triangle.wandzura_xiao_4(), 1.0e-14),
    (quadpy.triangle.wandzura_xiao_5(), 1.0e-14),
    (quadpy.triangle.wandzura_xiao_6(), 1.0e-14),
    (quadpy.triangle.williams_shunn_jameson_1(), 1.0e-14),
    (quadpy.triangle.williams_shunn_jameson_2(), 1.0e-14),
    (quadpy.triangle.williams_shunn_jameson_3(), 1.0e-13),
    (quadpy.triangle.williams_shunn_jameson_4(), 1.0e-14),
    (quadpy.triangle.williams_shunn_jameson_5(), 1.0e-13),
    (quadpy.triangle.williams_shunn_jameson_6(), 1.0e-11),
    (quadpy.triangle.williams_shunn_jameson_7(), 1.0e-12),
    (quadpy.triangle.williams_shunn_jameson_8(), 1.0e-11),
    (quadpy.triangle.witherden_vincent_01(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_02(), 1.0e-14),
    # (quadpy.triangle.witherden_vincent_03(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_04(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_05(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_06(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_07(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_08(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_09(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_10(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_11(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_12(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_13(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_14(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_15(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_16(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_17(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_18(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_19(), 1.0e-14),
    (quadpy.triangle.witherden_vincent_20(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_01(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_02(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_03(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_04(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_05(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_06(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_07(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_08(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_09(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_10(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_11(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_12(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_13(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_14(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_15(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_16(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_17(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_18(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_19(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_20(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_21(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_22(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_23(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_24(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_25(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_26(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_27(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_28(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_29(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_30(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_31(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_32(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_33(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_34(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_35(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_36(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_37(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_38(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_39(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_40(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_41(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_42(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_43(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_44(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_45(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_46(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_47(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_48(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_49(), 1.0e-14),
    (quadpy.triangle.xiao_gimbutas_50(), 1.0e-14),
    (quadpy.triangle.zhang_cui_liu_1(), 1.0e-14),
    (quadpy.triangle.zhang_cui_liu_2(), 1.0e-14),
    (quadpy.triangle.zhang_cui_liu_3(), 1.0e-14),
]


def _integrate_exact(f, triangle):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference triangle [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the triangle. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    xi = sympy.DeferredVector("xi")
    x_xi = (
        +triangle[0] * (1 - xi[0] - xi[1]) + triangle[1] * xi[0] + triangle[2] * xi[1]
    )
    abs_det_J = 2 * quadpy.triangle.volume(triangle)
    exact = sympy.integrate(
        sympy.integrate(abs_det_J * f(x_xi), (xi[1], 0, 1 - xi[0])), (xi[0], 0, 1)
    )
    return float(exact)


@pytest.mark.parametrize("scheme,tol", schemes_tol)
def test_scheme(scheme, tol):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    print(scheme.name)

    def eval_orthopolys(x):
        bary = numpy.array([x[0], x[1], 1.0 - x[0] - x[1]])
        out = numpy.concatenate(
            orthopy.triangle.tree(bary, scheme.degree + 1, "normal", symbolic=False)
        )
        return out

    vals = scheme.integrate(eval_orthopolys, triangle)
    # Put vals back into the tree structure:
    # len(approximate[k]) == k+1
    approximate = [
        vals[k * (k + 1) // 2 : (k + 1) * (k + 2) // 2]
        for k in range(scheme.degree + 2)
    ]

    exact = [numpy.zeros(k + 1) for k in range(scheme.degree + 2)]
    exact[0][0] = numpy.sqrt(2.0) / 2

    degree = check_degree_ortho(approximate, exact, abs_tol=tol)

    assert degree >= scheme.degree, "{} -- Observed: {}, expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.triangle.xiao_gimbutas_10()])
def test_show(scheme):
    triangle = numpy.array(
        [
            [numpy.cos(0.5 * numpy.pi), numpy.sin(0.5 * numpy.pi)],
            [numpy.cos(7.0 / 6.0 * numpy.pi), numpy.sin(7.0 / 6.0 * numpy.pi)],
            [numpy.cos(11.0 / 6.0 * numpy.pi), numpy.sin(11.0 / 6.0 * numpy.pi)],
        ]
    )
    scheme.show(triangle)
    return


def test_volume():
    # Assert computation of triangle volume in 3D is correct
    triangle = numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.7, 0.4, 1.1]])
    ref = numpy.sqrt(3.0) / 2.0
    assert abs(quadpy.triangle.get_vol(triangle) - ref) < 1.0e-14 * ref

    triangle = numpy.array([[0.0, 0.0, 0.0], [0.3, 0.4, 0.5], [0.7, 0.4, 1.1]])
    ref = numpy.sqrt(0.0209)
    assert abs(quadpy.triangle.get_vol(triangle) - ref) < 1.0e-14 * ref
    return


if __name__ == "__main__":
    # scheme_ = quadpy.triangle.WandzuraXiao(3)
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    from helpers import find_equal

    schemes_ = [scheme[0] for scheme in schemes_tol]
    find_equal(schemes_)
