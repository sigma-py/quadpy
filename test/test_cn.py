import numpy as np
import orthopy
import pytest

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.cn.cools_haegemans_1(n) for n in range(1, 7)]
    + [quadpy.cn.cools_haegemans_2(n) for n in range(2, 7)]
    + [quadpy.cn.dobrodeev_1970(n) for n in range(5, 8)]
    + [quadpy.cn.dobrodeev_1978(n) for n in range(2, 8)]
    + [quadpy.cn.hammer_stroud_1n(n) for n in range(2, 7)]
    + [quadpy.cn.hammer_stroud_2n(n) for n in range(2, 7)]
    + [quadpy.cn.mcnamee_stenger_3(n) for n in range(2, 7)]
    + [quadpy.cn.mcnamee_stenger_5(n) for n in range(2, 7)]
    + [quadpy.cn.mcnamee_stenger_7a(n) for n in range(3, 7)]
    + [quadpy.cn.mcnamee_stenger_7b(n) for n in range(3, 7)]
    + [quadpy.cn.mcnamee_stenger_9a(n) for n in range(4, 7)]
    + [quadpy.cn.mcnamee_stenger_9b(n) for n in range(4, 7)]
    + [quadpy.cn.stroud_cn_1_1(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_1_2(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_2_1(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_2_2(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_3_1(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_3_2(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_3_3(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_3_4(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_3_5(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_3_6(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_2(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_3(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_4(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_5(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_6(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_7(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_5_8(n) for n in range(3, 7)]
    + [quadpy.cn.stroud_cn_5_9(n) for n in range(2, 7)]
    + [quadpy.cn.stroud_cn_7_1(n) for n in range(2, 7)],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    print(scheme)

    n = scheme.dim
    cn_limits = [[-1.0, 1.0]] * n
    cn = quadpy.cn.ncube_points(*cn_limits)

    evaluator = orthopy.cn.Eval(scheme.points)

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), cn)
        exact = evaluator.int_p0 * 2 ** n if k == 0 else 0.0
        err = np.abs(approximate - exact)
        if np.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    max_err = np.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


# https://github.com/nschloe/quadpy/issues/401
def test_multiple_cubes_volume():
    val = quadpy.c2.schemes["stroud_c2_5_4"]().integrate(
        lambda x: x[0],
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
    assert val.shape == (5,)
    assert np.all(np.abs(val - 0.5) < 1.0e-13)


if __name__ == "__main__":
    test_multiple_cubes_volume()
