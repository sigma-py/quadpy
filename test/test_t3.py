import numpy
import pytest
import sympy

import quadpy
from helpers import check_degree
from quadpy.tn._helpers import integrate_monomial_over_unit_simplex

schemes = [
    quadpy.t3.beckers_haegemans_8(),
    quadpy.t3.beckers_haegemans_9(),
    quadpy.t3.gatermann(),
    quadpy.t3.hammer_marlowe_stroud_1(),
    quadpy.t3.hammer_marlowe_stroud_2(),
    quadpy.t3.hammer_marlowe_stroud_3(),
    quadpy.t3.hammer_stroud_2(),
    quadpy.t3.hammer_stroud_3(),
    quadpy.t3.jaskowiec_sukumar_02(),
    quadpy.t3.jaskowiec_sukumar_03(),
    quadpy.t3.jaskowiec_sukumar_04(),
    quadpy.t3.jaskowiec_sukumar_05(),
    quadpy.t3.jaskowiec_sukumar_06(),
    quadpy.t3.jaskowiec_sukumar_07(),
    quadpy.t3.jaskowiec_sukumar_08(),
    quadpy.t3.jaskowiec_sukumar_09(),
    quadpy.t3.jaskowiec_sukumar_10(),
    quadpy.t3.jaskowiec_sukumar_11(),
    quadpy.t3.jaskowiec_sukumar_12(),
    quadpy.t3.jaskowiec_sukumar_13(),
    quadpy.t3.jaskowiec_sukumar_14(),
    quadpy.t3.jaskowiec_sukumar_15(),
    quadpy.t3.jaskowiec_sukumar_16(),
    quadpy.t3.jaskowiec_sukumar_17(),
    quadpy.t3.jaskowiec_sukumar_18(),
    quadpy.t3.jaskowiec_sukumar_19a(),
    quadpy.t3.jaskowiec_sukumar_19b(),
    quadpy.t3.jaskowiec_sukumar_20(),
    quadpy.t3.keast_0(),
    quadpy.t3.keast_1(),
    quadpy.t3.keast_2(),
    quadpy.t3.keast_3(),
    quadpy.t3.keast_4(),
    quadpy.t3.keast_5(),
    quadpy.t3.keast_6(),
    quadpy.t3.keast_7(),
    quadpy.t3.keast_8(),
    quadpy.t3.keast_9(),
    quadpy.t3.liu_vinokur_01(),
    quadpy.t3.liu_vinokur_02(),
    quadpy.t3.liu_vinokur_03(),
    quadpy.t3.liu_vinokur_04(),
    quadpy.t3.liu_vinokur_05(),
    quadpy.t3.liu_vinokur_06(),
    quadpy.t3.liu_vinokur_07(),
    quadpy.t3.liu_vinokur_08(),
    quadpy.t3.liu_vinokur_09(),
    quadpy.t3.liu_vinokur_10(),
    quadpy.t3.liu_vinokur_11(),
    quadpy.t3.liu_vinokur_12(),
    quadpy.t3.liu_vinokur_13(),
    quadpy.t3.liu_vinokur_14(),
    quadpy.t3.maeztu_sainz(),
    quadpy.t3.stroud_t3_5_1(),
    quadpy.t3.stroud_t3_7_1(),
    quadpy.t3.shunn_ham_1(),
    quadpy.t3.shunn_ham_2(),
    quadpy.t3.shunn_ham_3(),
    quadpy.t3.shunn_ham_4(),
    quadpy.t3.shunn_ham_5(),
    quadpy.t3.shunn_ham_6(),
    quadpy.t3.vioreanu_rokhlin_0(),
    quadpy.t3.vioreanu_rokhlin_1(),
    quadpy.t3.vioreanu_rokhlin_2(),
    quadpy.t3.vioreanu_rokhlin_3(),
    quadpy.t3.vioreanu_rokhlin_4(),
    quadpy.t3.vioreanu_rokhlin_5(),
    quadpy.t3.vioreanu_rokhlin_6(),
    quadpy.t3.vioreanu_rokhlin_7(),
    quadpy.t3.vioreanu_rokhlin_8(),
    quadpy.t3.vioreanu_rokhlin_9(),
    quadpy.t3.xiao_gimbutas_01(),
    quadpy.t3.xiao_gimbutas_02(),
    quadpy.t3.xiao_gimbutas_03(),
    quadpy.t3.xiao_gimbutas_04(),
    quadpy.t3.xiao_gimbutas_05(),
    quadpy.t3.xiao_gimbutas_06(),
    quadpy.t3.xiao_gimbutas_07(),
    quadpy.t3.xiao_gimbutas_08(),
    quadpy.t3.xiao_gimbutas_09(),
    quadpy.t3.xiao_gimbutas_10(),
    quadpy.t3.xiao_gimbutas_11(),
    quadpy.t3.xiao_gimbutas_12(),
    quadpy.t3.xiao_gimbutas_13(),
    quadpy.t3.xiao_gimbutas_14(),
    quadpy.t3.xiao_gimbutas_15(),
    quadpy.t3.yu_1(),
    quadpy.t3.yu_2(),
    quadpy.t3.yu_3(),
    quadpy.t3.yu_4(),
    quadpy.t3.yu_5(),
    quadpy.t3.zhang_cui_liu_1(),
    quadpy.t3.zhang_cui_liu_2(),
    quadpy.t3.walkington_p5(),
    quadpy.t3.williams_shunn_jameson(),
    quadpy.t3.witherden_vincent_01(),
    quadpy.t3.witherden_vincent_02(),
    quadpy.t3.witherden_vincent_03(),
    quadpy.t3.witherden_vincent_05(),
    quadpy.t3.witherden_vincent_06(),
    quadpy.t3.witherden_vincent_07(),
    quadpy.t3.witherden_vincent_08(),
    quadpy.t3.witherden_vincent_09(),
    quadpy.t3.witherden_vincent_10(),
]


def _integrate_exact(f, t3):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference t3 [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the t3. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    xi = sympy.DeferredVector("xi")
    x_xi = (
        +t3[0] * (1 - xi[0] - xi[1] - xi[2])
        + t3[1] * xi[0]
        + t3[2] * xi[1]
        + t3[3] * xi[2]
    )
    abs_det_J = 6 * quadpy.t3.volume(t3)
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

    print(scheme)

    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    t3 = numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, t3),
        integrate_monomial_over_unit_simplex,
        3,
        scheme.degree + 1,
        scheme.test_tolerance,
    )

    assert (
        degree >= scheme.degree
    ), "{} -- observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.skip(reason="gh-actions's python cannot use system vtk")
@pytest.mark.parametrize("scheme", [quadpy.t3.hammer_marlowe_stroud_3()])
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


if __name__ == "__main__":
    # scheme_ = quadpy.t3.Stroud("T3 7-1")
    # test_scheme(scheme_)
    # # test_show(scheme_)
    # quadpy.t3.show(scheme_, backend="vtk")
    from helpers import find_equal

    find_equal(schemes)
