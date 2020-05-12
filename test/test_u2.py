import matplotlib.pyplot as plt
import numpy
import pytest

import quadpy
from helpers import check_degree
from quadpy.sn._helpers import integrate_monomial_over_unit_nsphere


@pytest.mark.parametrize("scheme", [quadpy.u2.krylov(k) for k in range(1, 6)])
def test_scheme(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: scheme.integrate(poly, [0.0, 0.0], 1.0),
        integrate_monomial_over_unit_nsphere,
        2,
        scheme.degree + 1,
    )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize("scheme", [quadpy.u2.krylov(3)])
def test_show(scheme):
    scheme.show(scheme)
    plt.close()
    return


if __name__ == "__main__":
    scheme_ = quadpy.u2.Krylov(30)
    test_scheme(scheme_)
    test_show(scheme_)
