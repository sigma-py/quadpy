import math

import numpy
import pytest

import orthopy
import quadpy
from helpers import check_degree_ortho


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

    # degree = check_degree(
    #     lambda poly: scheme.integrate(poly, cn),
    #     lambda exp: integrate_monomial_over_cn(cn_limits, exp),
    #     n,
    #     scheme.degree + 1,
    #     tol=tol,
    # )
    # assert degree >= scheme.degree, "observed: {}, expected: {}".format(
    #     degree, scheme.degree
    # )

    def eval_orthopolys(x):
        return numpy.concatenate(
            orthopy.ncube.tree(x, scheme.degree + 1, symbolic=False)
        )

    vals = scheme.integrate(eval_orthopolys, cn)

    # Put vals back into the tree structure:
    # len(approximate[k]) == k+1
    approximate = [
        vals[
            numpy.prod(range(k, k + n))
            // math.factorial(n) : numpy.prod(range(k + 1, k + 1 + n))
            // math.factorial(n)
        ]
        for k in range(scheme.degree + 2)
    ]

    exact = [numpy.zeros(len(s)) for s in approximate]
    exact[0][0] = numpy.sqrt(2.0) ** n

    degree, err = check_degree_ortho(approximate, exact, abs_tol=scheme.test_tolerance)

    assert (
        degree >= scheme.degree
    ), "{} (dim={})  --  observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, n, degree, scheme.degree, err
    )


if __name__ == "__main__":
    n_ = 4
    scheme_ = quadpy.cn.Stroud(n_, "Cn 7-1")
    test_scheme(scheme_, 1.0e-14)
