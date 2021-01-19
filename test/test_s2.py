import numpy as np
import orthopy
import pytest
from helpers import find_best_scheme

import quadpy


@pytest.mark.parametrize("scheme", quadpy.s2.schemes.values())
def test_scheme(scheme):
    try:
        scheme = scheme()
    except TypeError:
        scheme = scheme(1)

    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    evaluator = orthopy.s2.xu.Eval(scheme.points, "normal")

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), [0.0, 0.0], 1.0)
        exact = evaluator.int_p0 if k == 0 else 0.0
        err = np.abs(approximate - exact)
        max_err = max(max_err, np.max(err))
        if np.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), [0.0, 0.0], 1.0)
            exact = evaluator.int_p0 if i == 0 else 0.0
            err = np.abs(approximate - exact)
            max_err = max(max_err, np.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
        )


@pytest.mark.parametrize("scheme", [quadpy.s2.schemes["lether"](3)])
def test_show(scheme):
    scheme.show()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.s2.schemes.values(),
            degree,
            lambda pts: np.all(pts[0] ** 2 + pts[1] ** 2 <= 1),
            lambda keys: "plain" not in keys,
        )
        if best is None:
            break

        b = quadpy.s2.get_good_scheme(degree)

        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 20


if __name__ == "__main__":
    test_get_good_scheme()
