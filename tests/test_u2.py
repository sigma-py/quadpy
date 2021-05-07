import ndim
import numpy as np
import pytest
from helpers import check_degree
from matplotlib import pyplot as plt

import quadpy


@pytest.mark.parametrize(
    "scheme", [quadpy.u2.schemes["krylov"](k) for k in range(1, 6)]
)
def test_scheme(scheme):
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, [0.0, 0.0], 1.0),
        ndim.nsphere.integrate_monomial,
        2,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{}  --  observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.u2.schemes["krylov"](3)])
def test_show(scheme):
    scheme.show(scheme)
    plt.close()


if __name__ == "__main__":
    scheme_ = quadpy.u2.Krylov(30)
    test_scheme(scheme_)
    test_show(scheme_)
