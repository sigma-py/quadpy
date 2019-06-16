# -*- coding: utf-8 -*-
#
import numpy
import pytest
import accupy

import quadpy

from helpers import check_degree, integrate_monomial_over_enr


@pytest.mark.parametrize(
    "scheme,tol",
    [(scheme(n), 1.0e-14) for n in range(4, 6) for scheme in quadpy.enr.Stroud.values()]
    + [
        (scheme(n), 1.0e-14)
        for n in range(2, 6)
        for scheme in quadpy.enr.StroudSecrest.values()
    ],
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    n = scheme.dim
    degree = check_degree(
        lambda poly: quadpy.enr.integrate(poly, scheme, dot=accupy.fdot),
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
