import numpy
import orthopy
import pytest

import quadpy

schemes = [
    quadpy.t2.albrecht_collatz(),
    quadpy.t2.centroid(),
    quadpy.t2.cools_haegemans_1(),
    # quadpy.t2.cools_haegemans_2(),
    quadpy.t2.cubtri(),
    quadpy.t2.berntsen_espelid_1(),
    quadpy.t2.berntsen_espelid_2(),
    quadpy.t2.berntsen_espelid_3(),
    quadpy.t2.berntsen_espelid_4(),
    quadpy.t2.dcutri(),
    quadpy.t2.dunavant_01(),
    quadpy.t2.dunavant_02(),
    quadpy.t2.dunavant_03(),
    quadpy.t2.dunavant_04(),
    quadpy.t2.dunavant_05(),
    quadpy.t2.dunavant_06(),
    quadpy.t2.dunavant_07(),
    quadpy.t2.dunavant_08(),
    quadpy.t2.dunavant_09(),
    quadpy.t2.dunavant_10(),
    quadpy.t2.dunavant_11(),
    quadpy.t2.dunavant_12(),
    quadpy.t2.dunavant_13(),
    quadpy.t2.dunavant_14(),
    quadpy.t2.dunavant_15(),
    quadpy.t2.dunavant_16(),
    quadpy.t2.dunavant_17(),
    quadpy.t2.dunavant_18(),
    quadpy.t2.dunavant_19(),
    quadpy.t2.dunavant_20(),
    quadpy.t2.franke_09(),
    quadpy.t2.franke_10(),
    quadpy.t2.gatermann(),
    quadpy.t2.griener_schmid_1(),
    quadpy.t2.griener_schmid_2(),
    quadpy.t2.hammer_marlowe_stroud_1(),
    quadpy.t2.hammer_marlowe_stroud_2(),
    quadpy.t2.hammer_marlowe_stroud_3(),
    quadpy.t2.hammer_marlowe_stroud_4(),
    quadpy.t2.hammer_marlowe_stroud_5(),
    quadpy.t2.hillion_01(),
    quadpy.t2.hillion_02(),
    quadpy.t2.hillion_03(),
    quadpy.t2.hillion_04(),
    quadpy.t2.hillion_05(),
    quadpy.t2.hillion_06(),
    quadpy.t2.hillion_07(),
    quadpy.t2.hillion_08(),
    quadpy.t2.hillion_09(),
    quadpy.t2.hillion_10(),
    quadpy.t2.laursen_gellert_01(),
    quadpy.t2.laursen_gellert_02a(),
    quadpy.t2.laursen_gellert_02b(),
    quadpy.t2.laursen_gellert_03(),
    quadpy.t2.laursen_gellert_04(),
    quadpy.t2.laursen_gellert_05(),
    quadpy.t2.laursen_gellert_06(),
    quadpy.t2.laursen_gellert_07(),
    quadpy.t2.laursen_gellert_08(),
    quadpy.t2.laursen_gellert_09(),
    quadpy.t2.laursen_gellert_10(),
    quadpy.t2.laursen_gellert_11(),
    quadpy.t2.laursen_gellert_12(),
    quadpy.t2.laursen_gellert_13(),
    quadpy.t2.laursen_gellert_14(),
    quadpy.t2.laursen_gellert_15a(),
    quadpy.t2.laursen_gellert_15b(),
    quadpy.t2.lether(1),
    quadpy.t2.lether(2),
    quadpy.t2.lether(3),
    quadpy.t2.lether(4),
    quadpy.t2.lether(5),
    quadpy.t2.lether(6),
    quadpy.t2.lether(7),
    quadpy.t2.lether(8),
    quadpy.t2.lether(9),
    quadpy.t2.liu_vinokur_01(),
    quadpy.t2.liu_vinokur_02(),
    quadpy.t2.liu_vinokur_03(),
    quadpy.t2.liu_vinokur_04(),
    quadpy.t2.liu_vinokur_05(),
    quadpy.t2.liu_vinokur_06(),
    quadpy.t2.liu_vinokur_07(),
    quadpy.t2.liu_vinokur_08(),
    quadpy.t2.liu_vinokur_09(),
    quadpy.t2.liu_vinokur_10(),
    quadpy.t2.liu_vinokur_11(),
    quadpy.t2.liu_vinokur_12(),
    quadpy.t2.liu_vinokur_13(),
    quadpy.t2.lyness_jespersen_01(),
    quadpy.t2.lyness_jespersen_02(),
    quadpy.t2.lyness_jespersen_03(),
    quadpy.t2.lyness_jespersen_04(),
    quadpy.t2.lyness_jespersen_05(),
    quadpy.t2.lyness_jespersen_06(),
    quadpy.t2.lyness_jespersen_07(),
    quadpy.t2.lyness_jespersen_08(),
    quadpy.t2.lyness_jespersen_09(),
    quadpy.t2.lyness_jespersen_10(),
    quadpy.t2.lyness_jespersen_11(),
    quadpy.t2.lyness_jespersen_12(),
    quadpy.t2.lyness_jespersen_13(),
    quadpy.t2.lyness_jespersen_14(),
    quadpy.t2.lyness_jespersen_15(),
    quadpy.t2.lyness_jespersen_16(),
    quadpy.t2.lyness_jespersen_17(),
    quadpy.t2.lyness_jespersen_18(),
    quadpy.t2.lyness_jespersen_19(),
    quadpy.t2.lyness_jespersen_20(),
    quadpy.t2.lyness_jespersen_21(),
    quadpy.t2.papanicolopulos_sym_0(),
    quadpy.t2.papanicolopulos_sym_1(),
    quadpy.t2.papanicolopulos_sym_2(),
    quadpy.t2.papanicolopulos_sym_3(),
    quadpy.t2.papanicolopulos_sym_4(),
    quadpy.t2.papanicolopulos_sym_5(),
    quadpy.t2.papanicolopulos_sym_6(),
    quadpy.t2.papanicolopulos_sym_7(),
    quadpy.t2.papanicolopulos_sym_8(),
    quadpy.t2.papanicolopulos_rot_08(),
    quadpy.t2.papanicolopulos_rot_09(),
    quadpy.t2.papanicolopulos_rot_10(),
    quadpy.t2.papanicolopulos_rot_11(),
    quadpy.t2.papanicolopulos_rot_12(),
    quadpy.t2.papanicolopulos_rot_13(),
    quadpy.t2.papanicolopulos_rot_14(),
    quadpy.t2.papanicolopulos_rot_15(),
    quadpy.t2.papanicolopulos_rot_16(),
    quadpy.t2.papanicolopulos_rot_17(),
    quadpy.t2.seven_point(),
    quadpy.t2.strang_fix_cowper_01(),
    quadpy.t2.strang_fix_cowper_02(),
    quadpy.t2.strang_fix_cowper_03(),
    quadpy.t2.strang_fix_cowper_04(),
    quadpy.t2.strang_fix_cowper_05(),
    quadpy.t2.strang_fix_cowper_06(),
    quadpy.t2.strang_fix_cowper_07(),
    quadpy.t2.strang_fix_cowper_08(),
    quadpy.t2.strang_fix_cowper_09(),
    quadpy.t2.strang_fix_cowper_10(),
    quadpy.t2.stroud_t2_3_1(),
    quadpy.t2.stroud_t2_5_1(),
    quadpy.t2.stroud_t2_7_1(),
    quadpy.t2.taylor_wingate_bos_1(),
    quadpy.t2.taylor_wingate_bos_2(),
    quadpy.t2.taylor_wingate_bos_4(),
    quadpy.t2.taylor_wingate_bos_5(),
    quadpy.t2.taylor_wingate_bos_8(),
    quadpy.t2.vertex(),
    quadpy.t2.vioreanu_rokhlin_00(),
    quadpy.t2.vioreanu_rokhlin_01(),
    quadpy.t2.vioreanu_rokhlin_02(),
    quadpy.t2.vioreanu_rokhlin_03(),
    quadpy.t2.vioreanu_rokhlin_04(),
    quadpy.t2.vioreanu_rokhlin_05(),
    quadpy.t2.vioreanu_rokhlin_06(),
    quadpy.t2.vioreanu_rokhlin_07(),
    quadpy.t2.vioreanu_rokhlin_08(),
    quadpy.t2.vioreanu_rokhlin_09(),
    quadpy.t2.vioreanu_rokhlin_10(),
    quadpy.t2.vioreanu_rokhlin_11(),
    quadpy.t2.vioreanu_rokhlin_12(),
    quadpy.t2.vioreanu_rokhlin_13(),
    quadpy.t2.vioreanu_rokhlin_14(),
    quadpy.t2.vioreanu_rokhlin_15(),
    quadpy.t2.vioreanu_rokhlin_16(),
    quadpy.t2.vioreanu_rokhlin_17(),
    quadpy.t2.vioreanu_rokhlin_18(),
    quadpy.t2.vioreanu_rokhlin_19(),
    quadpy.t2.walkington_p5(),
    quadpy.t2.wandzura_xiao_1(),
    quadpy.t2.wandzura_xiao_2(),
    quadpy.t2.wandzura_xiao_3(),
    quadpy.t2.wandzura_xiao_4(),
    quadpy.t2.wandzura_xiao_5(),
    quadpy.t2.wandzura_xiao_6(),
    quadpy.t2.williams_shunn_jameson_1(),
    quadpy.t2.williams_shunn_jameson_2(),
    quadpy.t2.williams_shunn_jameson_3(),
    quadpy.t2.williams_shunn_jameson_4(),
    quadpy.t2.williams_shunn_jameson_5(),
    quadpy.t2.williams_shunn_jameson_6(),
    quadpy.t2.williams_shunn_jameson_7(),
    quadpy.t2.williams_shunn_jameson_8(),
    quadpy.t2.witherden_vincent_01(),
    quadpy.t2.witherden_vincent_02(),
    # quadpy.t2.witherden_vincent_03(),
    quadpy.t2.witherden_vincent_04(),
    quadpy.t2.witherden_vincent_05(),
    quadpy.t2.witherden_vincent_06(),
    quadpy.t2.witherden_vincent_07(),
    quadpy.t2.witherden_vincent_08(),
    quadpy.t2.witherden_vincent_09(),
    quadpy.t2.witherden_vincent_10(),
    quadpy.t2.witherden_vincent_11(),
    quadpy.t2.witherden_vincent_12(),
    quadpy.t2.witherden_vincent_13(),
    quadpy.t2.witherden_vincent_14(),
    quadpy.t2.witherden_vincent_15(),
    quadpy.t2.witherden_vincent_16(),
    quadpy.t2.witherden_vincent_17(),
    quadpy.t2.witherden_vincent_18(),
    quadpy.t2.witherden_vincent_19(),
    quadpy.t2.witherden_vincent_20(),
    quadpy.t2.xiao_gimbutas_01(),
    quadpy.t2.xiao_gimbutas_02(),
    quadpy.t2.xiao_gimbutas_03(),
    quadpy.t2.xiao_gimbutas_04(),
    quadpy.t2.xiao_gimbutas_05(),
    quadpy.t2.xiao_gimbutas_06(),
    quadpy.t2.xiao_gimbutas_07(),
    quadpy.t2.xiao_gimbutas_08(),
    quadpy.t2.xiao_gimbutas_09(),
    quadpy.t2.xiao_gimbutas_10(),
    quadpy.t2.xiao_gimbutas_11(),
    quadpy.t2.xiao_gimbutas_12(),
    quadpy.t2.xiao_gimbutas_13(),
    quadpy.t2.xiao_gimbutas_14(),
    quadpy.t2.xiao_gimbutas_15(),
    quadpy.t2.xiao_gimbutas_16(),
    quadpy.t2.xiao_gimbutas_17(),
    quadpy.t2.xiao_gimbutas_18(),
    quadpy.t2.xiao_gimbutas_19(),
    quadpy.t2.xiao_gimbutas_20(),
    quadpy.t2.xiao_gimbutas_21(),
    quadpy.t2.xiao_gimbutas_22(),
    quadpy.t2.xiao_gimbutas_23(),
    quadpy.t2.xiao_gimbutas_24(),
    quadpy.t2.xiao_gimbutas_25(),
    quadpy.t2.xiao_gimbutas_26(),
    quadpy.t2.xiao_gimbutas_27(),
    quadpy.t2.xiao_gimbutas_28(),
    quadpy.t2.xiao_gimbutas_29(),
    quadpy.t2.xiao_gimbutas_30(),
    quadpy.t2.xiao_gimbutas_31(),
    quadpy.t2.xiao_gimbutas_32(),
    quadpy.t2.xiao_gimbutas_33(),
    quadpy.t2.xiao_gimbutas_34(),
    quadpy.t2.xiao_gimbutas_35(),
    quadpy.t2.xiao_gimbutas_36(),
    quadpy.t2.xiao_gimbutas_37(),
    quadpy.t2.xiao_gimbutas_38(),
    quadpy.t2.xiao_gimbutas_39(),
    quadpy.t2.xiao_gimbutas_40(),
    quadpy.t2.xiao_gimbutas_41(),
    quadpy.t2.xiao_gimbutas_42(),
    quadpy.t2.xiao_gimbutas_43(),
    quadpy.t2.xiao_gimbutas_44(),
    quadpy.t2.xiao_gimbutas_45(),
    quadpy.t2.xiao_gimbutas_46(),
    quadpy.t2.xiao_gimbutas_47(),
    quadpy.t2.xiao_gimbutas_48(),
    quadpy.t2.xiao_gimbutas_49(),
    quadpy.t2.xiao_gimbutas_50(),
    quadpy.t2.zhang_cui_liu_1(),
    quadpy.t2.zhang_cui_liu_2(),
    quadpy.t2.zhang_cui_liu_3(),
]


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])

    evaluator = orthopy.t2.Eval(scheme.points, "normal")

    # assert contiguous x
    def f(x):
        assert x.flags["C_CONTIGUOUS"]
        assert x.shape[0] == 2
        return numpy.ones(x.shape[1:])

    approximate = scheme.integrate(f, triangle)

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), triangle)
        exact = evaluator.int_p0 * 0.5 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        max_err = max(max_err, numpy.max(err))
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), triangle)
            exact = evaluator.int_p0 * 0.5 if i == 0 else 0.0
            err = numpy.abs(approximate - exact)
            max_err = max(max_err, numpy.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
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


def test_get_good_scheme():
    for degree in range(51):
        best = None
        for scheme in schemes:
            if scheme.degree < degree:
                continue

            # allow only positive weights
            if any(scheme.weights < 0):
                continue

            # disallow points outside of the domain
            if numpy.any(scheme.points < 0):
                continue

            keys = set(scheme.symmetry_data.keys())

            if len(keys - set(["s1", "s2", "s3", "vertex"])) > 0:
                continue

            if best is not None:
                if len(scheme.weights) > len(best.weights):
                    continue
                if len(scheme.weights) == len(best.weights):
                    ratio = max(numpy.abs(scheme.weights)) / min(
                        numpy.abs(scheme.weights)
                    )
                    bratio = max(numpy.abs(best.weights)) / min(numpy.abs(best.weights))
                    if ratio > bratio:
                        continue
                    # check if it's actually the same scheme
                    if numpy.all(numpy.abs(scheme.points - best.points) < 1.0e-12):
                        continue

            # okay, looks like we found a better one!
            best = scheme

        print(degree, best.name)
        # print(best)
    return


if __name__ == "__main__":
    test_get_good_scheme()
    # test_multidim()
    # scheme_ = quadpy.t2.WandzuraXiao(3)
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    # from helpers import find_equal
    # schemes_ = [scheme[0] for scheme in schemes_tol]
    # find_equal(schemes_)
