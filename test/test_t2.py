import numpy
import pytest

# from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
import sympy

import orthopy
import quadpy
from helpers import check_degree_ortho

schemes_tol = [
    (quadpy.t2.albrecht_collatz(), 1.0e-14),
    (quadpy.t2.centroid(), 1.0e-14),
    (quadpy.t2.cools_haegemans_1(), 1.0e-13),
    # (quadpy.t2.cools_haegemans_2(), 1.0e-13),
    (quadpy.t2.cubtri(), 1.0e-14),
    (quadpy.t2.berntsen_espelid_1(), 1.0e-11),
    (quadpy.t2.berntsen_espelid_2(), 1.0e-11),
    (quadpy.t2.berntsen_espelid_3(), 1.0e-11),
    (quadpy.t2.berntsen_espelid_4(), 1.0e-11),
    (quadpy.t2.dcutri(), 1.0e-11),
    (quadpy.t2.dunavant_01(), 1.0e-14),
    (quadpy.t2.dunavant_02(), 1.0e-14),
    (quadpy.t2.dunavant_03(), 1.0e-14),
    (quadpy.t2.dunavant_04(), 1.0e-14),
    (quadpy.t2.dunavant_05(), 1.0e-12),
    (quadpy.t2.dunavant_06(), 1.0e-14),
    (quadpy.t2.dunavant_07(), 1.0e-14),
    (quadpy.t2.dunavant_08(), 1.0e-14),
    (quadpy.t2.dunavant_09(), 1.0e-14),
    (quadpy.t2.dunavant_10(), 1.0e-13),
    (quadpy.t2.dunavant_11(), 1.0e-12),
    (quadpy.t2.dunavant_12(), 1.0e-14),
    (quadpy.t2.dunavant_13(), 1.0e-14),
    (quadpy.t2.dunavant_14(), 1.0e-14),
    (quadpy.t2.dunavant_15(), 1.0e-13),
    (quadpy.t2.dunavant_16(), 1.0e-13),
    (quadpy.t2.dunavant_17(), 1.0e-13),
    (quadpy.t2.dunavant_18(), 1.0e-12),
    (quadpy.t2.dunavant_19(), 1.0e-13),
    (quadpy.t2.dunavant_20(), 1.0e-13),
    (quadpy.t2.franke_09(), 1.0e-14),
    (quadpy.t2.franke_10(), 1.0e-14),
    (quadpy.t2.gatermann(), 1.0e-12),
    (quadpy.t2.griener_schmid_1(), 1.0e-14),
    (quadpy.t2.griener_schmid_2(), 1.0e-14),
    (quadpy.t2.hammer_marlowe_stroud_1(), 1.0e-14),
    (quadpy.t2.hammer_marlowe_stroud_2(), 1.0e-14),
    (quadpy.t2.hammer_marlowe_stroud_3(), 1.0e-14),
    (quadpy.t2.hammer_marlowe_stroud_4(), 1.0e-14),
    (quadpy.t2.hammer_marlowe_stroud_5(), 1.0e-14),
    (quadpy.t2.hillion_01(), 1.0e-14),
    (quadpy.t2.hillion_02(), 1.0e-14),
    (quadpy.t2.hillion_03(), 1.0e-14),
    (quadpy.t2.hillion_04(), 1.0e-14),
    (quadpy.t2.hillion_05(), 1.0e-14),
    (quadpy.t2.hillion_06(), 1.0e-14),
    (quadpy.t2.hillion_07(), 1.0e-14),
    (quadpy.t2.hillion_08(), 1.0e-14),
    (quadpy.t2.hillion_09(), 1.0e-14),
    (quadpy.t2.hillion_10(), 1.0e-14),
    (quadpy.t2.laursen_gellert_01(), 1.0e-14),
    (quadpy.t2.laursen_gellert_02a(), 1.0e-14),
    (quadpy.t2.laursen_gellert_02b(), 1.0e-14),
    (quadpy.t2.laursen_gellert_03(), 1.0e-14),
    (quadpy.t2.laursen_gellert_04(), 1.0e-14),
    (quadpy.t2.laursen_gellert_05(), 1.0e-14),
    (quadpy.t2.laursen_gellert_06(), 1.0e-14),
    (quadpy.t2.laursen_gellert_07(), 1.0e-14),
    (quadpy.t2.laursen_gellert_08(), 1.0e-14),
    (quadpy.t2.laursen_gellert_09(), 1.0e-14),
    (quadpy.t2.laursen_gellert_10(), 1.0e-14),
    (quadpy.t2.laursen_gellert_11(), 1.0e-13),
    (quadpy.t2.laursen_gellert_12(), 1.0e-14),
    (quadpy.t2.laursen_gellert_13(), 1.0e-14),
    (quadpy.t2.laursen_gellert_14(), 1.0e-14),
    (quadpy.t2.laursen_gellert_15a(), 1.0e-14),
    (quadpy.t2.laursen_gellert_15b(), 1.0e-14),
    (quadpy.t2.lether(1), 1.0e-14),
    (quadpy.t2.lether(2), 1.0e-14),
    (quadpy.t2.lether(3), 1.0e-14),
    (quadpy.t2.lether(4), 1.0e-14),
    (quadpy.t2.lether(5), 1.0e-14),
    (quadpy.t2.lether(6), 1.0e-14),
    (quadpy.t2.lether(7), 1.0e-14),
    (quadpy.t2.lether(8), 1.0e-14),
    (quadpy.t2.lether(9), 1.0e-14),
    (quadpy.t2.liu_vinokur_01(), 1.0e-14),
    (quadpy.t2.liu_vinokur_02(), 1.0e-14),
    (quadpy.t2.liu_vinokur_03(), 1.0e-14),
    (quadpy.t2.liu_vinokur_04(), 1.0e-14),
    (quadpy.t2.liu_vinokur_05(), 1.0e-14),
    (quadpy.t2.liu_vinokur_06(), 1.0e-14),
    (quadpy.t2.liu_vinokur_07(), 1.0e-14),
    (quadpy.t2.liu_vinokur_08(), 1.0e-14),
    (quadpy.t2.liu_vinokur_09(), 1.0e-14),
    (quadpy.t2.liu_vinokur_10(), 1.0e-14),
    (quadpy.t2.liu_vinokur_11(), 1.0e-14),
    (quadpy.t2.liu_vinokur_12(), 1.0e-14),
    (quadpy.t2.liu_vinokur_13(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_01(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_02(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_03(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_04(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_05(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_06(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_07(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_08(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_09(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_10(), 1.0e-12),
    (quadpy.t2.lyness_jespersen_11(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_12(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_13(), 1.0e-13),
    (quadpy.t2.lyness_jespersen_14(), 1.0e-12),
    (quadpy.t2.lyness_jespersen_15(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_16(), 1.0e-11),
    (quadpy.t2.lyness_jespersen_17(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_18(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_19(), 1.0e-14),
    (quadpy.t2.lyness_jespersen_20(), 1.0e-12),
    (quadpy.t2.lyness_jespersen_21(), 1.0e-13),
    (quadpy.t2.papanicolopulos_sym_0(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_1(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_2(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_3(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_4(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_5(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_6(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_7(), 1.0e-14),
    (quadpy.t2.papanicolopulos_sym_8(), 1.0e-13),
    (quadpy.t2.papanicolopulos_rot_08(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_09(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_10(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_11(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_12(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_13(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_14(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_15(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_16(), 1.0e-14),
    (quadpy.t2.papanicolopulos_rot_17(), 1.0e-14),
    (quadpy.t2.seven_point(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_01(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_02(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_03(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_04(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_05(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_06(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_07(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_08(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_09(), 1.0e-14),
    (quadpy.t2.strang_fix_cowper_10(), 1.0e-14),
    (quadpy.t2.stroud_t2_3_1(), 1.0e-14),
    (quadpy.t2.stroud_t2_5_1(), 1.0e-14),
    (quadpy.t2.stroud_t2_7_1(), 1.0e-12),
    (quadpy.t2.taylor_wingate_bos_1(), 1.0e-14),
    (quadpy.t2.taylor_wingate_bos_2(), 1.0e-12),
    (quadpy.t2.taylor_wingate_bos_4(), 1.0e-12),
    (quadpy.t2.taylor_wingate_bos_5(), 1.0e-12),
    (quadpy.t2.taylor_wingate_bos_8(), 1.0e-12),
    (quadpy.t2.vertex(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_00(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_01(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_02(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_03(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_04(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_05(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_06(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_07(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_08(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_09(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_10(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_11(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_12(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_13(), 1.0e-13),
    (quadpy.t2.vioreanu_rokhlin_14(), 1.0e-13),
    (quadpy.t2.vioreanu_rokhlin_15(), 1.0e-13),
    (quadpy.t2.vioreanu_rokhlin_16(), 1.0e-12),
    (quadpy.t2.vioreanu_rokhlin_17(), 1.0e-14),
    (quadpy.t2.vioreanu_rokhlin_18(), 1.0e-13),
    (quadpy.t2.vioreanu_rokhlin_19(), 1.0e-11),
    (quadpy.t2.walkington_p5(), 1.0e-14),
    (quadpy.t2.wandzura_xiao_1(), 1.0e-14),
    (quadpy.t2.wandzura_xiao_2(), 1.0e-14),
    (quadpy.t2.wandzura_xiao_3(), 1.0e-14),
    (quadpy.t2.wandzura_xiao_4(), 1.0e-14),
    (quadpy.t2.wandzura_xiao_5(), 1.0e-14),
    (quadpy.t2.wandzura_xiao_6(), 1.0e-14),
    (quadpy.t2.williams_shunn_jameson_1(), 1.0e-14),
    (quadpy.t2.williams_shunn_jameson_2(), 1.0e-14),
    (quadpy.t2.williams_shunn_jameson_3(), 1.0e-13),
    (quadpy.t2.williams_shunn_jameson_4(), 1.0e-14),
    (quadpy.t2.williams_shunn_jameson_5(), 1.0e-13),
    (quadpy.t2.williams_shunn_jameson_6(), 1.0e-11),
    (quadpy.t2.williams_shunn_jameson_7(), 1.0e-12),
    (quadpy.t2.williams_shunn_jameson_8(), 1.0e-11),
    (quadpy.t2.witherden_vincent_01(), 1.0e-14),
    (quadpy.t2.witherden_vincent_02(), 1.0e-14),
    # (quadpy.t2.witherden_vincent_03(), 1.0e-14),
    (quadpy.t2.witherden_vincent_04(), 1.0e-14),
    (quadpy.t2.witherden_vincent_05(), 1.0e-14),
    (quadpy.t2.witherden_vincent_06(), 1.0e-14),
    (quadpy.t2.witherden_vincent_07(), 1.0e-14),
    (quadpy.t2.witherden_vincent_08(), 1.0e-14),
    (quadpy.t2.witherden_vincent_09(), 1.0e-14),
    (quadpy.t2.witherden_vincent_10(), 1.0e-14),
    (quadpy.t2.witherden_vincent_11(), 1.0e-14),
    (quadpy.t2.witherden_vincent_12(), 1.0e-14),
    (quadpy.t2.witherden_vincent_13(), 1.0e-14),
    (quadpy.t2.witherden_vincent_14(), 1.0e-14),
    (quadpy.t2.witherden_vincent_15(), 1.0e-14),
    (quadpy.t2.witherden_vincent_16(), 1.0e-14),
    (quadpy.t2.witherden_vincent_17(), 1.0e-14),
    (quadpy.t2.witherden_vincent_18(), 1.0e-14),
    (quadpy.t2.witherden_vincent_19(), 1.0e-14),
    (quadpy.t2.witherden_vincent_20(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_01(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_02(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_03(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_04(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_05(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_06(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_07(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_08(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_09(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_10(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_11(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_12(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_13(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_14(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_15(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_16(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_17(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_18(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_19(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_20(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_21(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_22(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_23(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_24(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_25(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_26(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_27(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_28(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_29(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_30(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_31(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_32(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_33(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_34(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_35(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_36(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_37(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_38(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_39(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_40(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_41(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_42(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_43(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_44(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_45(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_46(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_47(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_48(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_49(), 1.0e-14),
    (quadpy.t2.xiao_gimbutas_50(), 1.0e-14),
    (quadpy.t2.zhang_cui_liu_1(), 1.0e-14),
    (quadpy.t2.zhang_cui_liu_2(), 1.0e-14),
    (quadpy.t2.zhang_cui_liu_3(), 1.0e-14),
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
    abs_det_J = 2 * quadpy.t2.volume(triangle)
    exact = sympy.integrate(
        sympy.integrate(abs_det_J * f(x_xi), (xi[1], 0, 1 - xi[0])), (xi[0], 0, 1)
    )
    return float(exact)


@pytest.mark.parametrize("scheme,tol", schemes_tol)
def test_scheme(scheme, tol):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])

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


@pytest.mark.parametrize("scheme", [quadpy.t2.xiao_gimbutas_10()])
def test_show(scheme):
    triangle = numpy.array(
        [
            [numpy.cos(0.5 * numpy.pi), numpy.sin(0.5 * numpy.pi)],
            [numpy.cos(7.0 / 6.0 * numpy.pi), numpy.sin(7.0 / 6.0 * numpy.pi)],
            [numpy.cos(11.0 / 6.0 * numpy.pi), numpy.sin(11.0 / 6.0 * numpy.pi)],
        ]
    )
    scheme.show(triangle)


def test_volume():
    # Assert computation of triangle volume in 3D is correct
    triangle = numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.7, 0.4, 1.1]])
    ref = numpy.sqrt(3.0) / 2.0
    assert abs(quadpy.t2.get_vol(triangle) - ref) < 1.0e-14 * ref

    triangle = numpy.array([[0.0, 0.0, 0.0], [0.3, 0.4, 0.5], [0.7, 0.4, 1.1]])
    ref = numpy.sqrt(0.0209)
    assert abs(quadpy.t2.get_vol(triangle) - ref) < 1.0e-14 * ref


def test_multidim():
    scheme = quadpy.t2.dunavant_05()

    numpy.random.seed(0)
    # simple scalar integration
    tri = numpy.random.rand(3, 2)
    val = scheme.integrate(lambda x: numpy.sin(x[0]), tri)
    assert val.shape == ()

    # scalar integration on 4 subdomains
    tri = numpy.random.rand(3, 4, 2)
    val = scheme.integrate(lambda x: numpy.sin(x[0]), tri)
    assert val.shape == (4,)

    # scalar integration in 4D
    tri = numpy.random.rand(3, 4)
    val = scheme.integrate(lambda x: numpy.sin(x[0]), tri)
    assert val.shape == ()

    # vector-valued integration on 4 subdomains
    tri = numpy.random.rand(3, 4, 2)
    val = scheme.integrate(lambda x: [numpy.sin(x[0]), numpy.cos(x[1])], tri)
    assert val.shape == (2, 4)

    # vector-valued integration in 4D
    tri = numpy.random.rand(3, 4)
    val = scheme.integrate(lambda x: [numpy.sin(x[0]), numpy.cos(x[1])], tri)
    assert val.shape == (2,)

    # # another vector-valued integration in 3D
    # # This is one case where the integration routine may not properly recognize the
    # # dimensionality of the domain. Use the `dim` parameter.
    # val = scheme.integrate(
    #     lambda x: [
    #         x[0] + numpy.sin(x[1]),
    #         numpy.cos(x[0]) * x[2],
    #         numpy.sin(x[0]) + x[1] + x[2],
    #     ],
    #     [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
    #     dim=1,
    # )
    # assert val.shape == (3,)


if __name__ == "__main__":
    test_multidim()
    # scheme_ = quadpy.t2.WandzuraXiao(3)
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    # from helpers import find_equal
    # schemes_ = [scheme[0] for scheme in schemes_tol]
    # find_equal(schemes_)
