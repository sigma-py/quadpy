import numpy as np
import orthopy
import pytest
from helpers import find_best_scheme
from matplotlib import pyplot as plt

import quadpy


@pytest.mark.parametrize("scheme", quadpy.e3r2.schemes.values())
def test_scheme(scheme, tol=1.0e-14):
    scheme = scheme()

    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    evaluator = orthopy.enr2.Eval(scheme.points, "physicists")

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator))
        exact = evaluator.int_p0 if k == 0 else 0.0
        err = np.abs(approximate - exact)
        if np.any(err > tol):
            break
        k += 1

    max_err = np.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e3r2.schemes["stroud_secrest_07"]()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.e3r2.schemes.values(),
            degree,
            lambda pts: True,
            lambda keys: len(
                keys - {"zero3", "symm_r00", "symm_rr0", "symm_rrr", "symm_rrs"}
            )
            == 0,
        )
        if best is None:
            break

        # print(degree, best.name)
        b = quadpy.e3r2.get_good_scheme(degree)
        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 8


if __name__ == "__main__":
    test_get_good_scheme()
