import numpy as np

import quadpy


def test_u2():
    scheme = quadpy.u2.schemes["krylov"](3)
    scheme.integrate(lambda x: np.exp(x[0]), np.array([0.0, 0.3]), 0.7)
    scheme = quadpy.u2.schemes["krylov"](5)
    scheme.integrate(
        lambda x: [np.exp(x[0]), np.exp(x[0])],
        np.array([[1.0, 1.0], [0.0, 0.3], [2.0, 2.0]]),
        [1.0, 0.7, 0.333],
    )


def test_s2():
    scheme = quadpy.s2.get_good_scheme(5)
    scheme.integrate(lambda x: np.exp(x[0]), np.array([0.0, 0.3]), 0.7)
    scheme = quadpy.s2.get_good_scheme(6)
    scheme.integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.array([[1.0, 1.0], [0.0, 0.3], [2.0, 2.0]]),
        [1.0, 0.7, 0.333],
    )


def test_c3():
    scheme = quadpy.c3.product(quadpy.c1.newton_cotes_closed(3))
    val = scheme.integrate(
        lambda x: np.exp(x[0]),
        quadpy.c3.cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
    )

    scheme = quadpy.c3.product(quadpy.c1.newton_cotes_closed(3))
    val = scheme.integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.stack(
            [
                quadpy.c3.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.c3.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.c3.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.c3.cube_points([0, 1], [0, 1], [0, 1]),
                quadpy.c3.cube_points([0, 1], [0, 1], [0, 1]),
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 5)


def test_c1():
    quadpy.c1.newton_cotes_closed(4).integrate(np.exp, [0.0, 1.0])
    quadpy.c1.newton_cotes_closed(4).integrate(np.exp, np.array([[0.0], [1.0]]))
    quadpy.c1.newton_cotes_closed(4).integrate(
        lambda x: [np.exp(x), np.sin(x), np.cos(x)],
        np.array([[0.0, 1.0], [1.0, 2.0]]),
    )


def test_p3():
    scheme = quadpy.p3.felippa_3()
    scheme.integrate(
        lambda x: np.exp(x[0]),
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]]),
    )

    scheme = quadpy.p3.felippa_5()
    scheme.integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.array(
            [
                [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]],
                [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
                [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]],
                [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
            ]
        ),
    )


def test_c2():
    quadpy.c2.schemes["stroud_c2_5_4"]().integrate(
        lambda x: np.exp(x[0]),
        quadpy.c2.rectangle_points([0.0, 1.0], [0.0, 1.0]),
    )

    val = quadpy.c2.schemes["stroud_c2_3_1"]().integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.stack(
            [
                quadpy.c2.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.c2.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.c2.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.c2.rectangle_points([0.0, 1.0], [0.0, 1.0]),
                quadpy.c2.rectangle_points([0.0, 1.0], [0.0, 1.0]),
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 5)


def test_e2r2():
    scheme = quadpy.e2r2.get_good_scheme(3)
    scheme.integrate(lambda x: np.exp(x[0]))


def test_e2r():
    scheme = quadpy.e2r.get_good_scheme(5)
    scheme.integrate(lambda x: np.exp(x[0]))


def test_u3():
    quadpy.u3.schemes["lebedev_003a"]().integrate(
        lambda x: np.exp(1j * x[0]) + 1j * x[0] ** 2,
        np.array([0.0, 0.3, 0.0]),
        0.7,
    )

    # quadpy.u3.lebedev_003a().integrate(
    #     lambda x: [np.exp(x[0]), np.exp(x[1])],
    #     np.array([[1.0, 1.0, 0.0], [0.0, 0.3, 0.0], [2.0, 2.0, 0.0]]),
    #     [1.0, 0.7, 0.333],
    # )


def test_s3():
    scheme = quadpy.s3.get_good_scheme(3)
    scheme.integrate(lambda x: np.exp(x[0]), [0.0, 0.0, 0.0], 1.0)

    scheme = quadpy.s3.get_good_scheme(5)
    scheme.integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        [[1.0, 1.0, 0.0], [0.0, 0.3, 0.0], [2.0, 2.0, 0.0]],
        [1.0, 0.7, 0.333],
    )


def test_t3():
    scheme = quadpy.t3.get_good_scheme(4)
    scheme.integrate(
        lambda x: np.exp(x[0]),
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
    )

    scheme.integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.stack(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0], [0.0, 0, 1]],
            ],
            axis=-2,
        ),
    )


def test_t2():
    quadpy.t2.schemes["cubtri"]().integrate(
        lambda x: np.exp(x[0]), np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    )

    val = quadpy.t2.schemes["cubtri"]().integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.stack(
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


def test_w3():
    quadpy.w3.felippa_4().integrate(
        lambda x: np.exp(x[0]),
        np.array(
            [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]]
        ),
    )

    val = quadpy.w3.felippa_4().integrate(
        lambda x: [np.exp(x[0]), np.exp(x[1])],
        np.stack(
            [
                [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
                [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
                [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
            ],
            axis=-2,
        ),
    )
    assert val.shape == (2, 3)


def test_tn():
    dim = 4
    quadpy.tn.grundmann_moeller(dim, 3).integrate(
        lambda x: np.exp(x[0]),
        np.array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [1.0, 2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 3.0, 1.0, 0.0],
                [0.0, 0.0, 4.0, 1.0],
            ]
        ),
    )


def test_sn():
    dim = 4
    quadpy.sn.dobrodeev_1970(dim).integrate(lambda x: np.exp(x[0]), np.zeros(4), 1.0)


def test_cn():
    dim = 4
    quadpy.cn.stroud_cn_3_3(dim).integrate(
        lambda x: np.exp(x[0]),
        quadpy.cn.ncube_points([0.0, 1.0], [0.1, 0.9], [-1.0, 1.0], [-1.0, -0.5]),
    )


def test_un():
    dim = 4
    quadpy.un.dobrodeev_1978(dim).integrate(lambda x: np.exp(x[0]), np.zeros(dim), 1.0)


def test_enr2():
    dim = 4
    quadpy.enr2.stroud_enr2_5_2(dim).integrate(lambda x: np.exp(x[0]))


def test_e1r():
    scheme = quadpy.e1r.gauss_laguerre(5)
    scheme.integrate(lambda x: x[0] ** 2)


def test_e3r():
    quadpy.e3r.get_good_scheme(5).integrate(lambda x: np.exp(x[0]))


def test_e1r2():
    scheme = quadpy.e1r2.gauss_hermite(5)
    scheme.integrate(lambda x: x[0] ** 2)


if __name__ == "__main__":
    test_u3()
