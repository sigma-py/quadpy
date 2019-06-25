# -*- coding: utf-8 -*-
#
import numpy
import pytest
import accupy

import quadpy

from helpers import check_degree, integrate_monomial_over_enr


@pytest.mark.parametrize(
    "scheme",
    [quadpy.enr.stroud_enr_3_1(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_3_2(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_5_1(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_5_3(n) for n in range(4, 6)]
    + [quadpy.enr.stroud_enr_5_4(n) for n in range(4, 6)]
    #
    + [quadpy.enr.stroud_secrest_i(n) for n in range(2, 6)]
    + [quadpy.enr.stroud_secrest_ii(n) for n in range(2, 6)]
    + [quadpy.enr.stroud_secrest_ii(n) for n in range(2, 6)]
    + [quadpy.enr.stroud_secrest_iv(n) for n in range(2, 6)],
)
def test_scheme(scheme, tol=1.0e-14):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    n = scheme.dim
    degree = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        integrate_monomial_over_enr,
        n,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


if __name__ == "__main__":
    dim_ = 2
    # quadpy.e3r2.show(quadpy.enr.Stroud(dim_, '5-1a'), backend='vtk')
    scheme_ = quadpy.enr.Stroud(dim_, "5-3")
    test_scheme(scheme_, 1.0e-14)
