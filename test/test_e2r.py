import accupy
import ndim
import numpy as np
import pytest
from helpers import check_degree, find_best_scheme

import quadpy


@pytest.mark.parametrize("scheme", quadpy.e2r.schemes.values())
def test_scheme(scheme):
    scheme = scheme()
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        ndim.enr.integrate_monomial,
        2,
        scheme.degree + 1,
        tol=scheme.test_tolerance * 1.1,
    )
    assert degree >= scheme.degree, (
        f"{scheme.name} -- observed: {degree}, expected: {scheme.degree} "
        f"(max err: {err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e2r.schemes["rabinowitz_richter_1"]()])
def test_show(scheme):
    scheme.show()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.e2r.schemes.values(),
            degree,
            lambda pts: True,
            lambda keys: "plain" not in keys,
        )
        if best is None:
            break

        b = quadpy.e2r.get_good_scheme(degree)

        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 16


if __name__ == "__main__":
    test_get_good_scheme()
