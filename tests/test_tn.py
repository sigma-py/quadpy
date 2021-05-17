import ndim
import numpy as np
import pytest
from helpers import check_degree

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.tn.grundmann_moeller(dim, k) for dim in range(3, 7) for k in range(5)]
    #
    + [
        quadpy.tn.silvester(dim, variant, n)
        for dim in range(1, 6)
        for variant in ["open", "closed"]
        for n in range(1, 7)
    ]
    #
    + [quadpy.tn.stroud_tn_1_1(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_1_2(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_2_1a(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_2_1b(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_2_2(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_1(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_2(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_3(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_4(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_5(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_6a(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_6b(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_7(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_8(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_9(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_4_1(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_5_1(dim) for dim in range(3, 7)]
    + [quadpy.tn.stroud_tn_3_10(dim) for dim in [3, 4, 6, 7]]
    + [quadpy.tn.stroud_tn_3_11(dim) for dim in [3, 4, 6, 7]]
    + [quadpy.tn.stroud_tn_5_2(dim) for dim in range(4, 8)]
    #
    + [quadpy.tn.walkington_1(dim) for dim in range(2, 8)]
    + [quadpy.tn.walkington_2(dim) for dim in range(2, 8)]
    + [quadpy.tn.walkington_3(dim) for dim in range(2, 8)]
    + [quadpy.tn.walkington_5(dim) for dim in [2, 3]]
    + [quadpy.tn.walkington_7(dim) for dim in [3]],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    print(scheme)

    n = scheme.dim
    simplex = np.zeros((n + 1, n))
    for k in range(n):
        simplex[k + 1, k] = 1
    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, simplex),
        ndim.nsimplex.integrate_monomial,
        n,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{} -- observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


if __name__ == "__main__":
    scheme_ = quadpy.tn.silvester(3, "open", 5)
    test_scheme(scheme_)
