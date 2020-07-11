import numpy
import orthopy
import pytest

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.cn.cools_haegemans_1(n) for n in range(1, 7)]
    + [quadpy.cn.cools_haegemans_2(n) for n in range(2, 7)]
    # [quadpy.cn.dobrodeev_1970(n) for n in range(5, 8)]
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
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    n = scheme.dim
    cn_limits = [[-1.0, 1.0]] * n
    cn = quadpy.cn.ncube_points(*cn_limits)

    evaluator = orthopy.cn.Eval(scheme.points.T)

    degree = None
    for k in range(scheme.degree + 2):
        approximate = scheme.integrate(lambda x: next(evaluator)[0], cn)
        exact = numpy.sqrt(2.0) ** n if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        if numpy.any(err > scheme.test_tolerance):
            degree = k - 1
            break

    max_err = numpy.max(err)
    assert degree >= scheme.degree, (
        f"{scheme.name} -- observed: {degree}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


if __name__ == "__main__":
    n_ = 4
    scheme_ = quadpy.cn.Stroud(n_, "Cn 7-1")
    test_scheme(scheme_, 1.0e-14)
