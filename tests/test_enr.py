import accupy
import ndim
import numpy as np
import pytest
from helpers import check_degree

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.enr.mcnamee_stenger_3(n) for n in range(1, 6)]
    + [quadpy.enr.mcnamee_stenger_5(n) for n in range(2, 6)]
    + [quadpy.enr.mcnamee_stenger_7a(n) for n in range(3, 6)]
    + [quadpy.enr.mcnamee_stenger_7b(n) for n in range(3, 6)]
    # The condition of the degree-9 schemes is so bad that the tolerence had to be
    # 1.0e-2
    # + [quadpy.enr.mcnamee_stenger_9a(n) for n in range(4, 6)]
    # + [quadpy.enr.mcnamee_stenger_9b(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_3_1(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_3_2(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_5_1(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_5_3(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_5_4(n) for n in range(4, 6)]
    #
    + [quadpy.enr.stroud_secrest_1(n) for n in range(2, 6)]
    + [quadpy.enr.stroud_secrest_2(n) for n in range(2, 6)]
    + [quadpy.enr.stroud_secrest_3(n) for n in range(2, 6)]
    + [quadpy.enr.stroud_secrest_4(n) for n in range(2, 6)],
)
def test_scheme(scheme):
    assert scheme.points.dtype == np.float64, scheme.name
    assert scheme.weights.dtype == np.float64, scheme.name

    print(scheme)

    n = scheme.dim
    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        ndim.enr.integrate_monomial,
        n,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{} (dim={}) -- Observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, n, degree, scheme.degree, err
    )


if __name__ == "__main__":
    dim_ = 2
    # quadpy.e3r2.show(quadpy.enr.Stroud(dim_, '5-1a'), backend='vtk')
    scheme_ = quadpy.enr.Stroud(dim_, "5-3")
    test_scheme(scheme_, 1.0e-14)
