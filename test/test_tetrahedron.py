import matplotlib.pyplot as plt
import numpy
import pytest
import sympy

import quadpy
from helpers import check_degree
from quadpy.nsimplex._helpers import integrate_monomial_over_unit_simplex

schemes = [
    quadpy.tetrahedron.beckers_haegemans_8(),
    quadpy.tetrahedron.beckers_haegemans_9(),
    quadpy.tetrahedron.gatermann(),
    quadpy.tetrahedron.hammer_marlowe_stroud_1(),
    quadpy.tetrahedron.hammer_marlowe_stroud_2(),
    quadpy.tetrahedron.hammer_marlowe_stroud_3(),
    quadpy.tetrahedron.hammer_stroud_2(),
    quadpy.tetrahedron.hammer_stroud_3(),
    quadpy.tetrahedron.keast_0(),
    quadpy.tetrahedron.keast_1(),
    quadpy.tetrahedron.keast_2(),
    quadpy.tetrahedron.keast_3(),
    quadpy.tetrahedron.keast_4(),
    quadpy.tetrahedron.keast_5(),
    quadpy.tetrahedron.keast_6(),
    quadpy.tetrahedron.keast_7(),
    quadpy.tetrahedron.keast_8(),
    quadpy.tetrahedron.keast_9(),
    quadpy.tetrahedron.liu_vinokur_01(),
    quadpy.tetrahedron.liu_vinokur_02(),
    quadpy.tetrahedron.liu_vinokur_03(),
    quadpy.tetrahedron.liu_vinokur_04(),
    quadpy.tetrahedron.liu_vinokur_05(),
    quadpy.tetrahedron.liu_vinokur_06(),
    quadpy.tetrahedron.liu_vinokur_07(),
    quadpy.tetrahedron.liu_vinokur_08(),
    quadpy.tetrahedron.liu_vinokur_09(),
    quadpy.tetrahedron.liu_vinokur_10(),
    quadpy.tetrahedron.liu_vinokur_11(),
    quadpy.tetrahedron.liu_vinokur_12(),
    quadpy.tetrahedron.liu_vinokur_13(),
    quadpy.tetrahedron.liu_vinokur_14(),
    quadpy.tetrahedron.maeztu_sainz(),
    quadpy.tetrahedron.newton_cotes_closed(1),
    quadpy.tetrahedron.newton_cotes_closed(2),
    quadpy.tetrahedron.newton_cotes_closed(3),
    quadpy.tetrahedron.newton_cotes_closed(4),
    quadpy.tetrahedron.newton_cotes_closed(5),
    quadpy.tetrahedron.newton_cotes_closed(6),
    quadpy.tetrahedron.newton_cotes_open(1),
    quadpy.tetrahedron.newton_cotes_open(2),
    quadpy.tetrahedron.newton_cotes_open(3),
    quadpy.tetrahedron.newton_cotes_open(4),
    quadpy.tetrahedron.newton_cotes_open(5),
    quadpy.tetrahedron.newton_cotes_open(6),
    quadpy.tetrahedron.stroud_t3_5_1(),
    quadpy.tetrahedron.stroud_t3_7_1(),
    quadpy.tetrahedron.shunn_ham_1(),
    quadpy.tetrahedron.shunn_ham_2(),
    quadpy.tetrahedron.shunn_ham_3(),
    quadpy.tetrahedron.shunn_ham_4(),
    quadpy.tetrahedron.shunn_ham_5(),
    quadpy.tetrahedron.shunn_ham_6(),
    quadpy.tetrahedron.vioreanu_rokhlin_0(),
    quadpy.tetrahedron.vioreanu_rokhlin_1(),
    quadpy.tetrahedron.vioreanu_rokhlin_2(),
    quadpy.tetrahedron.vioreanu_rokhlin_3(),
    quadpy.tetrahedron.vioreanu_rokhlin_4(),
    quadpy.tetrahedron.vioreanu_rokhlin_5(),
    quadpy.tetrahedron.vioreanu_rokhlin_6(),
    quadpy.tetrahedron.vioreanu_rokhlin_7(),
    quadpy.tetrahedron.vioreanu_rokhlin_8(),
    quadpy.tetrahedron.vioreanu_rokhlin_9(),
    quadpy.tetrahedron.xiao_gimbutas_01(),
    quadpy.tetrahedron.xiao_gimbutas_02(),
    quadpy.tetrahedron.xiao_gimbutas_03(),
    quadpy.tetrahedron.xiao_gimbutas_04(),
    quadpy.tetrahedron.xiao_gimbutas_05(),
    quadpy.tetrahedron.xiao_gimbutas_06(),
    quadpy.tetrahedron.xiao_gimbutas_07(),
    quadpy.tetrahedron.xiao_gimbutas_08(),
    quadpy.tetrahedron.xiao_gimbutas_09(),
    quadpy.tetrahedron.xiao_gimbutas_10(),
    quadpy.tetrahedron.xiao_gimbutas_11(),
    quadpy.tetrahedron.xiao_gimbutas_12(),
    quadpy.tetrahedron.xiao_gimbutas_13(),
    quadpy.tetrahedron.xiao_gimbutas_14(),
    quadpy.tetrahedron.xiao_gimbutas_15(),
    quadpy.tetrahedron.yu_1(),
    quadpy.tetrahedron.yu_2(),
    quadpy.tetrahedron.yu_3(),
    quadpy.tetrahedron.yu_4(),
    quadpy.tetrahedron.yu_5(),
    quadpy.tetrahedron.zhang_cui_liu_1(),
    quadpy.tetrahedron.zhang_cui_liu_2(),
    quadpy.tetrahedron.walkington_p5(),
    quadpy.tetrahedron.williams_shunn_jameson(),
    quadpy.tetrahedron.witherden_vincent_01(),
    quadpy.tetrahedron.witherden_vincent_02(),
    quadpy.tetrahedron.witherden_vincent_03(),
    quadpy.tetrahedron.witherden_vincent_05(),
    quadpy.tetrahedron.witherden_vincent_06(),
    quadpy.tetrahedron.witherden_vincent_07(),
    quadpy.tetrahedron.witherden_vincent_08(),
    quadpy.tetrahedron.witherden_vincent_09(),
    quadpy.tetrahedron.witherden_vincent_10(),
]


def _integrate_exact(f, tetrahedron):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference tetrahedron [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the tetrahedron. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    xi = sympy.DeferredVector("xi")
    x_xi = (
        +tetrahedron[0] * (1 - xi[0] - xi[1] - xi[2])
        + tetrahedron[1] * xi[0]
        + tetrahedron[2] * xi[1]
        + tetrahedron[3] * xi[2]
    )
    abs_det_J = 6 * quadpy.tetrahedron.volume(tetrahedron)
    exact = sympy.integrate(
        sympy.integrate(
            sympy.integrate(abs_det_J * f(x_xi), (xi[2], 0, 1 - xi[0] - xi[1])),
            (xi[1], 0, 1 - xi[0]),
        ),
        (xi[0], 0, 1),
    )
    return float(exact)


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    tetrahedron = numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    degree = check_degree(
        lambda poly: scheme.integrate(poly, tetrahedron),
        integrate_monomial_over_unit_simplex,
        3,
        scheme.degree + 1,
    )
    assert degree == scheme.degree, "{} -- Observed: {}, expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.tetrahedron.hammer_marlowe_stroud_3()])
def test_show(scheme):
    tet = numpy.array(
        [
            [numpy.cos(0.5 * numpy.pi), numpy.sin(0.5 * numpy.pi), -0.5],
            [numpy.cos(7.0 / 6.0 * numpy.pi), numpy.sin(7.0 / 6.0 * numpy.pi), -0.5],
            [numpy.cos(11.0 / 6.0 * numpy.pi), numpy.sin(11.0 / 6.0 * numpy.pi), -0.5],
            [0.0, 0.0, 1.0],
        ]
    )
    scheme.show(tet, render=False)
    plt.close()
    return


if __name__ == "__main__":
    # scheme_ = quadpy.tetrahedron.Stroud("T3 7-1")
    # test_scheme(scheme_)
    # # test_show(scheme_)
    # quadpy.tetrahedron.show(scheme_, backend="vtk")
    from helpers import find_equal

    find_equal(schemes)
