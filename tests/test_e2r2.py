import numpy as np
import orthopy
import pytest
from helpers import find_best_scheme

import quadpy


@pytest.mark.parametrize("scheme", quadpy.e2r2.schemes.values())
def test_scheme(scheme):
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
        if np.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    max_err = np.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e2r2.schemes["rabinowitz_richter_1"]()])
def test_show(scheme):
    scheme.show()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.e2r2.schemes.values(),
            degree,
            lambda pts: True,
            lambda keys: "plain" not in keys,
        )
        if best is None:
            break

        b = quadpy.e2r2.get_good_scheme(degree)
        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 16


if __name__ == "__main__":
    test_get_good_scheme()
