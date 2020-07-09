import numpy
import orthopy
import pytest
from helpers import check_degree_ortho

import quadpy

schemes = [
    quadpy.s2.albrecht_1(),
    quadpy.s2.albrecht_2(),
    quadpy.s2.albrecht_3(),
    quadpy.s2.albrecht_4(),
    quadpy.s2.albrecht_5(),
    quadpy.s2.albrecht_6(),
    quadpy.s2.albrecht_7(),
    quadpy.s2.albrecht_8(),
    quadpy.s2.albrecht_collatz(),
    quadpy.s2.cools_haegemans_1(),
    quadpy.s2.cools_haegemans_2(),
    quadpy.s2.cools_haegemans_3(),
    quadpy.s2.cools_kim_1(),
    quadpy.s2.cools_kim_2(),
    quadpy.s2.cools_kim_3(),
    quadpy.s2.haegemans_piessens(),
    quadpy.s2.hammer_stroud_11_2(),
    quadpy.s2.hammer_stroud_12_2(),
    quadpy.s2.hammer_stroud_13_2(),
    quadpy.s2.hammer_stroud_17(),
    quadpy.s2.hammer_stroud_18(),
    quadpy.s2.hammer_stroud_19(),
    quadpy.s2.hammer_stroud_20(),
    quadpy.s2.hammer_stroud_21(),
    quadpy.s2.kim_song_1(),
    quadpy.s2.kim_song_2(),
    quadpy.s2.kim_song_3(),
    quadpy.s2.kim_song_4(),
    quadpy.s2.kim_song_5(),
    quadpy.s2.kim_song_6(),
    quadpy.s2.kim_song_7(),
    quadpy.s2.kim_song_8(),
    quadpy.s2.kim_song_9(),
    quadpy.s2.kim_song_10(),
    quadpy.s2.kim_song_11(),
    quadpy.s2.kim_song_12(),
    quadpy.s2.kim_song_13(),
    quadpy.s2.kim_song_14(),
    quadpy.s2.kim_song_15(),
    quadpy.s2.lether(1),
    quadpy.s2.lether(2),
    quadpy.s2.lether(3),
    quadpy.s2.lether(4),
    quadpy.s2.lether(5),
    quadpy.s2.lether(6),
    quadpy.s2.luo_meng_1(),
    quadpy.s2.luo_meng_2(),
    quadpy.s2.luo_meng_3(),
    quadpy.s2.luo_meng_4(),
    quadpy.s2.luo_meng_5(),
    # quadpy.s2.luo_meng_6(),
    # quadpy.s2.luo_meng_7(),
    quadpy.s2.mysovskih_1(),
    quadpy.s2.mysovskih_2(),
    quadpy.s2.mysovskih_3(),
    quadpy.s2.peirce_1956_1(),
    quadpy.s2.peirce_1956_2(),
    quadpy.s2.peirce_1956_3(),
    quadpy.s2.peirce_1957(1),
    quadpy.s2.peirce_1957(2),
    quadpy.s2.peirce_1957(3),
    quadpy.s2.peirce_1957(5),
    quadpy.s2.piessens_haegemans(),
    quadpy.s2.rabinowitz_richter_1(),
    quadpy.s2.rabinowitz_richter_2(),
    quadpy.s2.rabinowitz_richter_3(),
    quadpy.s2.rabinowitz_richter_4(),
    quadpy.s2.rabinowitz_richter_5(),
    quadpy.s2.rabinowitz_richter_6(),
    quadpy.s2.stroud_s2_3_1(),
    quadpy.s2.stroud_s2_3_2(),
    quadpy.s2.stroud_s2_4_1(),
    quadpy.s2.stroud_s2_5_1(),
    quadpy.s2.stroud_s2_5_2(),
    quadpy.s2.stroud_s2_7_1(),
    quadpy.s2.stroud_s2_7_2(),
    quadpy.s2.stroud_s2_9_1(),
    quadpy.s2.stroud_s2_9_2(),
    quadpy.s2.stroud_s2_9_3(),
    quadpy.s2.stroud_s2_9_4(),
    quadpy.s2.stroud_s2_9_5(),
    quadpy.s2.stroud_s2_11_1(),
    quadpy.s2.stroud_s2_11_2(),
    quadpy.s2.stroud_s2_11_3(),
    quadpy.s2.stroud_s2_11_4(),
    quadpy.s2.stroud_s2_13_1(),
    quadpy.s2.stroud_s2_13_2(),
    quadpy.s2.stroud_s2_15_1(),
    quadpy.s2.stroud_s2_15_2(),
    quadpy.s2.stroud_s2_17_1(),
    quadpy.s2.wissmann_becker_6_1(),
    quadpy.s2.wissmann_becker_6_2(),
    quadpy.s2.wissmann_becker_8_1(),
]


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    # degree = check_degree(
    #     lambda poly: scheme.integrate(poly, [0.0, 0.0], 1.0),
    #     integrate_monomial_over_unit_nball,
    #     2,
    #     scheme.degree + 1,
    #     tol=tol,
    # )
    # assert degree == scheme.degree, "{}  -- Observed: {}   expected: {}".format(
    #     scheme.name, degree, scheme.degree
    # )

    def eval_orthopolys(x):
        return numpy.concatenate(
            orthopy.disk.tree(x, scheme.degree + 1, symbolic=False)
        )

    vals = scheme.integrate(eval_orthopolys, [0, 0], 1)
    # Put vals back into the tree structure:
    # len(approximate[k]) == k+1
    approximate = [
        vals[k * (k + 1) // 2 : (k + 1) * (k + 2) // 2]
        for k in range(scheme.degree + 2)
    ]

    exact = [numpy.zeros(k + 1) for k in range(scheme.degree + 2)]
    exact[0][0] = numpy.sqrt(numpy.pi)

    degree, err = check_degree_ortho(approximate, exact, abs_tol=scheme.test_tolerance)

    assert (
        degree >= scheme.degree
    ), "{} -- Observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.s2.lether(3)])
def test_show(scheme):
    scheme.show()


if __name__ == "__main__":
    # scheme_ = quadpy.s2.lether(5)
    # scheme_ = quadpy.s2.albrecht_8()
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    from helpers import find_equal

    find_equal(schemes)
