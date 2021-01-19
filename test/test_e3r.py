import accupy
import ndim
import numpy as np
import pytest
from helpers import check_degree, find_best_scheme
from matplotlib import pyplot as plt

import quadpy


@pytest.mark.parametrize("scheme", quadpy.e3r.schemes.values())
def test_scheme(scheme):
    scheme = scheme()

    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        ndim.enr.integrate_monomial,
        3,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{} -- observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.e3r.schemes["stroud_secrest_10"]()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.e3r.schemes.values(),
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
        b = quadpy.e3r.get_good_scheme(degree)
        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 8


if __name__ == "__main__":
    scheme_ = quadpy.e3r.Stroud("5-1")
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_, backend="vtk")
