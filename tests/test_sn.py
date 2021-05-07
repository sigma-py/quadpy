import ndim
import numpy as np
import pytest
from helpers import check_degree

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.sn.dobrodeev_1970(n) for n in range(3, 9)]
    + [quadpy.sn.dobrodeev_1978(n) for n in range(2, 7)]
    + [quadpy.sn.mcnamee_stenger_3(n) for n in range(2, 7)]
    + [quadpy.sn.mcnamee_stenger_5(n) for n in range(2, 7)]
    + [quadpy.sn.mcnamee_stenger_7a(n) for n in range(3, 7)]
    + [quadpy.sn.mcnamee_stenger_7b(n) for n in range(3, 7)]
    + [quadpy.sn.mcnamee_stenger_9a(n) for n in range(4, 7)]
    + [quadpy.sn.mcnamee_stenger_9b(n) for n in range(4, 7)]
    + [quadpy.sn.stoyanova(n, variant_v_plus=True) for n in range(5, 10)]
    + [quadpy.sn.stoyanova(n, variant_v_plus=False) for n in range(5, 10)]
    + [quadpy.sn.stroud_sn_2_1(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_3_1(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_3_2(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_5_2(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_5_3(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_5_4(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_5_5(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_5_6(dim) for dim in range(2, 9)]
    + [quadpy.sn.stroud_sn_5_1a(dim) for dim in range(4, 8)]
    + [quadpy.sn.stroud_sn_5_1b(dim) for dim in range(4, 8)]
    + [quadpy.sn.stroud_sn_7_1a(dim) for dim in range(3, 8)]
    + [quadpy.sn.stroud_sn_7_1b(dim) for dim in range(3, 7)]
    + [quadpy.sn.stroud_sn_7_2(dim) for dim in range(3, 7)]
    + [quadpy.sn.stroud_sn_7_3a(dim) for dim in range(3, 7)]
    + [quadpy.sn.stroud_sn_7_3b(dim) for dim in range(3, 7)]
    + [quadpy.sn.stroud_sn_9_1a(dim) for dim in range(3, 7)]
    + [quadpy.sn.stroud_sn_9_1b(dim) for dim in range(4, 7)]
    + [quadpy.sn.stroud_sn_11_1a(dim) for dim in [3, 4]]
    + [quadpy.sn.stroud_sn_11_1b(dim) for dim in [4, 5]],
)
def test_scheme(scheme):
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    n = scheme.dim
    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, center=np.zeros(n), radius=1),
        ndim.nball.integrate_monomial,
        n,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert degree >= scheme.degree, (
        f"{scheme.name} -- observed: {degree}, expected: {scheme.degree} "
        f"(max err: {err:.3e})"
    )


if __name__ == "__main__":
    n_ = 3
    scheme_ = quadpy.sn.Stroud(n_, "Sn 2-1", symbolic=True)
    test_scheme(scheme_)
