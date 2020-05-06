import numpy
import pytest

import orthopy
import quadpy
from helpers import check_degree_ortho

schemes = [
    quadpy.disk.albrecht_1(),
    quadpy.disk.albrecht_2(),
    quadpy.disk.albrecht_3(),
    quadpy.disk.albrecht_4(),
    quadpy.disk.albrecht_5(),
    quadpy.disk.albrecht_6(),
    quadpy.disk.albrecht_7(),
    quadpy.disk.albrecht_8(),
    quadpy.disk.albrecht_collatz(),
    quadpy.disk.cools_haegemans_1(),
    quadpy.disk.cools_haegemans_2(),
    quadpy.disk.cools_haegemans_3(),
    quadpy.disk.cools_kim_1(),
    quadpy.disk.cools_kim_2(),
    quadpy.disk.cools_kim_3(),
    quadpy.disk.haegemans_piessens(),
    quadpy.disk.hammer_stroud_11_2(),
    quadpy.disk.hammer_stroud_12_2(),
    quadpy.disk.hammer_stroud_13_2(),
    quadpy.disk.hammer_stroud_17(),
    quadpy.disk.hammer_stroud_18(),
    quadpy.disk.hammer_stroud_19(),
    quadpy.disk.hammer_stroud_20(),
    quadpy.disk.hammer_stroud_21(),
    quadpy.disk.kim_song_1(),
    quadpy.disk.kim_song_2(),
    quadpy.disk.kim_song_3(),
    quadpy.disk.kim_song_4(),
    quadpy.disk.kim_song_5(),
    quadpy.disk.kim_song_6(),
    quadpy.disk.kim_song_7(),
    quadpy.disk.kim_song_8(),
    quadpy.disk.kim_song_9(),
    quadpy.disk.kim_song_10(),
    quadpy.disk.kim_song_11(),
    quadpy.disk.kim_song_12(),
    quadpy.disk.kim_song_13(),
    quadpy.disk.kim_song_14(),
    quadpy.disk.kim_song_15(),
    quadpy.disk.lether(1),
    quadpy.disk.lether(2),
    quadpy.disk.lether(3),
    quadpy.disk.lether(4),
    quadpy.disk.lether(5),
    quadpy.disk.lether(6),
    quadpy.disk.mysovskih_1(),
    quadpy.disk.mysovskih_2(),
    quadpy.disk.mysovskih_3(),
    quadpy.disk.peirce_1956_1(),
    quadpy.disk.peirce_1956_2(),
    quadpy.disk.peirce_1956_3(),
    quadpy.disk.peirce_1957(1),
    quadpy.disk.peirce_1957(2),
    quadpy.disk.peirce_1957(3),
    quadpy.disk.peirce_1957(5),
    quadpy.disk.piessens_haegemans(),
    quadpy.disk.rabinowitz_richter_1(),
    quadpy.disk.rabinowitz_richter_2(),
    quadpy.disk.rabinowitz_richter_3(),
    quadpy.disk.rabinowitz_richter_4(),
    quadpy.disk.rabinowitz_richter_5(),
    quadpy.disk.rabinowitz_richter_6(),
    quadpy.disk.stroud_s2_3_1(),
    quadpy.disk.stroud_s2_3_2(),
    quadpy.disk.stroud_s2_4_1(),
    quadpy.disk.stroud_s2_5_1(),
    quadpy.disk.stroud_s2_5_2(),
    quadpy.disk.stroud_s2_7_1(),
    quadpy.disk.stroud_s2_7_2(),
    quadpy.disk.stroud_s2_9_1(),
    quadpy.disk.stroud_s2_9_2(),
    quadpy.disk.stroud_s2_9_3(),
    quadpy.disk.stroud_s2_9_4(),
    quadpy.disk.stroud_s2_9_5(),
    quadpy.disk.stroud_s2_11_1(),
    quadpy.disk.stroud_s2_11_2(),
    quadpy.disk.stroud_s2_11_3(),
    quadpy.disk.stroud_s2_11_4(),
    quadpy.disk.stroud_s2_13_1(),
    quadpy.disk.stroud_s2_13_2(),
    quadpy.disk.stroud_s2_15_1(),
    quadpy.disk.stroud_s2_15_2(),
    quadpy.disk.stroud_s2_17_1(),
    quadpy.disk.wissmann_becker_6_1(),
    quadpy.disk.wissmann_becker_6_2(),
    quadpy.disk.wissmann_becker_8_1(),
]


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme, tol=1.0e-11):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

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

    degree = check_degree_ortho(approximate, exact, abs_tol=tol)

    assert degree >= scheme.degree, "{} -- Observed: {}, expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.disk.lether(3)])
def test_show(scheme):
    scheme.show()
    return


if __name__ == "__main__":
    # scheme_ = quadpy.disk.lether(5)
    # scheme_ = quadpy.disk.albrecht_8()
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    from helpers import find_equal

    find_equal(schemes)
