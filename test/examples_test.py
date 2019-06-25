import numpy

import quadpy


def test_circle():
    scheme = quadpy.circle.krylov(3)
    scheme.integrate(lambda x: numpy.exp(x[0]), numpy.array([0.0, 0.3]), 0.7)
    scheme = quadpy.circle.krylov(5)
    scheme.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[0])],
        numpy.array([[1.0, 1.0], [0.0, 0.3], [2.0, 2.0]]),
        [1.0, 0.7, 0.333],
    )
    return


def test_disk():
    scheme = quadpy.disk.peirce_1957(5)
    scheme.integrate(lambda x: numpy.exp(x[0]), numpy.array([0.0, 0.3]), 0.7)
    scheme = quadpy.disk.peirce_1957(5)
    scheme.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.array([[1.0, 1.0], [0.0, 0.3], [2.0, 2.0]]),
        [1.0, 0.7, 0.333],
    )
    return


def test_hexahedron():
    scheme = quadpy.hexahedron.product(quadpy.line_segment.newton_cotes_closed(3))
    val = scheme.integrate(
        lambda x: numpy.exp(x[0]),
        quadpy.hexahedron.cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
    )

    scheme = quadpy.hexahedron.product(quadpy.line_segment.newton_cotes_closed(3))
    val = scheme.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.stack(
            [
                quadpy.hexahedron.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.hexahedron.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.hexahedron.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.hexahedron.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.hexahedron.cube_points([0, 1], [0, 1], [0, 1]),
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 5)
    return


def test_line_segment():
    quadpy.line_segment.newton_cotes_closed(4).integrate(numpy.exp, [0.0, 1.0])
    quadpy.line_segment.newton_cotes_closed(4).integrate(
        numpy.exp, numpy.array([[0.0], [1.0]])
    )
    quadpy.line_segment.newton_cotes_closed(4).integrate(
        lambda x: [numpy.exp(x), numpy.sin(x), numpy.cos(x)],
        numpy.array([[0.0, 1.0], [1.0, 2.0]]),
    )
    return


def test_pyramid():
    scheme = quadpy.pyramid.felippa_3()
    scheme.integrate(
        lambda x: numpy.exp(x[0]),
        numpy.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]]),
    )

    scheme = quadpy.pyramid.felippa_5()
    scheme.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.array(
            [
                [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]],
                [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
                [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]],
                [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
            ]
        ),
    )
    return


def test_quadrilateral():
    quadpy.quadrilateral.stroud_c2_5_4().integrate(
        lambda x: numpy.exp(x[0]),
        quadpy.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0]),
    )

    val = quadpy.quadrilateral.stroud_c2_3_1().integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.stack(
            [
                quadpy.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.quadrilateral.rectangle_points([0.0, 1.0], [0.0, 1.0]),
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 5)
    return


def test_e2r2():
    scheme = quadpy.e2r2.rabinowitz_richter_3()
    scheme.integrate(lambda x: numpy.exp(x[0]))
    return


def test_e2r():
    scheme = quadpy.e2r.rabinowitz_richter_5()
    scheme.integrate(lambda x: numpy.exp(x[0]))
    return


def test_sphere():
    quadpy.sphere.lebedev_003a().integrate(
        lambda x: numpy.exp(1j * x[0]) + 1j * x[0] ** 2,
        numpy.array([0.0, 0.3, 0.0]),
        0.7,
    )

    quadpy.sphere.lebedev_003a().integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.array([[1.0, 1.0, 0.0], [0.0, 0.3, 0.0], [2.0, 2.0, 0.0]]),
        [1.0, 0.7, 0.333],
    )
    return


def test_ball():
    scheme = quadpy.ball.hammer_stroud_15_3a()
    scheme.integrate(lambda x: numpy.exp(x[0]), [0.0, 0.0, 0.0], 1.0)

    scheme = quadpy.ball.hammer_stroud_15_3b()
    scheme.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        [[1.0, 1.0, 0.0], [0.0, 0.3, 0.0], [2.0, 2.0, 0.0]],
        [1.0, 0.7, 0.333],
    )
    return


def test_tetrahedron():
    quadpy.tetrahedron.shunn_ham_3().integrate(
        lambda x: numpy.exp(x[0]),
        numpy.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
    )

    quadpy.tetrahedron.shunn_ham_3().integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.stack(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
            ],
            axis=-2,
        ),
    )
    return


def test_triangle():
    quadpy.triangle.cubtri().integrate(
        lambda x: numpy.exp(x[0]), numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    )

    val = quadpy.triangle.cubtri().integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.stack(
            [
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 5)
    return


def test_wedge():
    quadpy.wedge.felippa_4().integrate(
        lambda x: numpy.exp(x[0]),
        numpy.array(
            [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]]
        ),
    )

    val = quadpy.wedge.felippa_4().integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.stack(
            [
                [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
                [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
                [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 3)
    return


def test_nsimplex():
    dim = 4
    quadpy.nsimplex.grundmann_moeller(dim, 3).integrate(
        lambda x: numpy.exp(x[0]),
        numpy.array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [1.0, 2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 3.0, 1.0, 0.0],
                [0.0, 0.0, 4.0, 1.0],
            ]
        ),
    )
    return


def test_nball():
    dim = 4
    quadpy.nball.dobrodeev_1970(dim).integrate(
        lambda x: numpy.exp(x[0]), numpy.zeros(4), 1.0
    )
    return


def test_ncube():
    dim = 4
    quadpy.ncube.stroud_cn_3_3(dim).integrate(
        lambda x: numpy.exp(x[0]),
        quadpy.ncube.ncube_points([0.0, 1.0], [0.1, 0.9], [-1.0, 1.0], [-1.0, -0.5]),
    )
    return


def test_nsphere():
    dim = 4
    quadpy.nsphere.dobrodeev_1978(dim).integrate(
        lambda x: numpy.exp(x[0]), numpy.zeros(dim), 1.0
    )
    return


def test_enr2():
    dim = 4
    quadpy.enr2.stroud_enr2_5_2(dim).integrate(lambda x: numpy.exp(x[0]))
    return


def test_e1r():
    scheme = quadpy.e1r.gauss_laguerre(5)
    scheme.integrate(lambda x: x[0] ** 2)
    return


def test_e3r():
    quadpy.e3r.stroud_secrest_ix().integrate(lambda x: numpy.exp(x[0]))
    return


def test_e1r2():
    scheme = quadpy.e1r2.gauss_hermite(5)
    scheme.integrate(lambda x: x[0] ** 2)
    return


if __name__ == "__main__":
    test_sphere()
