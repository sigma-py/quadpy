import ndim
import numpy as np
import pytest
from helpers import check_degree

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.un.dobrodeev_1978(n) for n in range(2, 7)]
    + [quadpy.un.mysovskikh_1(n) for n in range(2, 7)]
    + [quadpy.un.mysovskikh_2(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_3_1(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_3_2(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_5_1(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_5_2(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_5_3(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_5_4(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_5_4(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_7_1(n) for n in range(2, 7)]
    + [quadpy.un.stroud_un_7_2(n) for n in range(2, 7)]
    # The scheme has degree 11, so don't push it too far with n. First of all,
    # the number of points increases exponentially, and so does the number of
    # polynomials of degree at most 11.
    + [quadpy.un.stroud_un_11_1(n) for n in range(3, 6)]
    + [quadpy.un.stroud_1967(n) for n in range(2, 7)],
)
def test_scheme(scheme):
    assert scheme.points.dtype in [np.int64, np.float64], scheme.name
    assert scheme.weights.dtype in [np.int64, np.float64], scheme.name

    print(scheme)

    n = scheme.dim
    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, center=np.zeros(n), radius=1),
        ndim.nsphere.integrate_monomial,
        n,
        scheme.degree + 1,
        tol=scheme.test_tolerance * 1.1,
    )
    assert degree >= scheme.degree, (
        f"{scheme.name} (dim={n})  --  "
        f"observed: {degree}, expected: {scheme.degree} (max err: {err:.3e})"
    )


if __name__ == "__main__":
    n_ = 5
    scheme_ = quadpy.un.Stroud(n_, "Un 11-1")
    test_scheme(scheme_)
