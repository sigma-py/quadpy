import numpy as np
import orthopy
import pytest

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.enr2.cools_haegemans_1(n) for n in range(1, 6)]
    + [quadpy.enr2.cools_haegemans_2(n) for n in range(1, 6)]
    + [quadpy.enr2.cools_haegemans_3(n) for n in range(2, 6)]
    + [quadpy.enr2.lu_darmofal_1(n) for n in range(2, 6)]
    + [quadpy.enr2.mcnamee_stenger_3(n) for n in range(1, 6)]
    + [quadpy.enr2.mcnamee_stenger_5(n) for n in range(2, 6)]
    + [quadpy.enr2.mcnamee_stenger_7a(n) for n in range(3, 6)]
    + [quadpy.enr2.mcnamee_stenger_7b(n) for n in range(3, 6)]
    + [quadpy.enr2.mcnamee_stenger_9a(n) for n in range(4, 6)]
    + [quadpy.enr2.mcnamee_stenger_9b(n) for n in range(4, 6)]
    + [quadpy.enr2.phillips(n) for n in range(2, 6)]
    + [quadpy.enr2.stenger_7a(n) for n in [3, 4]]
    + [quadpy.enr2.stenger_7b(n) for n in range(3, 6)]
    + [quadpy.enr2.stenger_9a(n) for n in range(3, 6)]
    + [quadpy.enr2.stenger_9b(n) for n in range(4, 6)]
    + [quadpy.enr2.stenger_11a(n) for n in range(3, 5)]
    + [quadpy.enr2.stenger_11b(n) for n in range(3, 6)]
    + [quadpy.enr2.stroud_enr2_3_1(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_enr2_3_2(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_enr2_5_1a(n) for n in range(2, 7)]
    + [quadpy.enr2.stroud_enr2_5_1b(n) for n in [3, 5, 6]]
    + [quadpy.enr2.stroud_enr2_5_2(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_enr2_5_3(n) for n in range(3, 8)]
    + [quadpy.enr2.stroud_enr2_5_4(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_enr2_5_5a(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_enr2_5_5b(n) for n in [2]]
    + [quadpy.enr2.stroud_enr2_5_6(n) for n in range(5, 8)]
    + [quadpy.enr2.stroud_enr2_7_1a(n) for n in [2, 3, 4, 6, 7]]
    + [quadpy.enr2.stroud_enr2_7_1b(n) for n in [3, 4]]
    + [quadpy.enr2.stroud_enr2_7_2(n) for n in range(3, 8)]
    + [quadpy.enr2.stroud_enr2_7_3a(n) for n in [3, 4]]
    + [quadpy.enr2.stroud_enr2_7_3b(n) for n in range(3, 6)]
    + [quadpy.enr2.stroud_enr2_9_1a(n) for n in range(3, 7)]
    + [quadpy.enr2.stroud_enr2_9_1b(n) for n in range(4, 7)]
    + [quadpy.enr2.stroud_enr2_11_1a(n) for n in range(3, 5)]
    + [quadpy.enr2.stroud_enr2_11_1b(n) for n in range(3, 6)]
    + [quadpy.enr2.stroud_1967_5_a(n) for n in range(2, 7)]
    + [quadpy.enr2.stroud_1967_5_b(n) for n in [3, 5, 6, 7]]
    + [quadpy.enr2.stroud_1967_7_2a(n) for n in [2, 3, 4, 6, 7]]
    + [quadpy.enr2.stroud_1967_7_2b(n) for n in [3, 4]]
    + [quadpy.enr2.stroud_1967_7_4(n) for n in range(3, 8)]
    + [quadpy.enr2.stroud_secrest_1(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_secrest_2(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_secrest_3(n) for n in range(2, 8)]
    + [quadpy.enr2.stroud_secrest_4(n) for n in range(2, 8)]
    + [quadpy.enr2.xiu(n) for n in range(2, 8)],
)
def test_scheme(scheme):
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    evaluator = orthopy.enr2.Eval(scheme.points, "physicists")

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator))
        exact = evaluator.int_p0 if k == 0 else 0.0
        err = np.abs(approximate - exact)
        if np.any(err > 1.0e-11):
            break
        k += 1

    max_err = np.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


if __name__ == "__main__":
    dim_ = 5
    # quadpy.e3r2.show(quadpy.enr2.Stroud(dim_, '5-1a'), backend='vtk')
    scheme_ = quadpy.enr2.Stroud(dim_, "11-1b")
    test_scheme(scheme_, 1.0e-14)
