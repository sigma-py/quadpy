# -*- coding: utf-8 -*-
#
import numpy
import pytest
import accupy

import quadpy

from helpers import check_degree, integrate_monomial_over_enr


@pytest.mark.parametrize(
    "scheme,tol",
    [(scheme(), 1.0e-14) for scheme in quadpy.e2r.HaegemansPiessens.values()]
    + [(scheme(), 1.0e-14) for scheme in quadpy.e2r.RabinowitzRichter.values()]
    + [(scheme(), 1.0e-14) for scheme in quadpy.e2r.Stroud.values()]
    + [(scheme(), 1.0e-14) for scheme in quadpy.e2r.StroudSecrest.values()],
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: quadpy.e2r.integrate(poly, scheme, dot=accupy.fdot),
        integrate_monomial_over_enr,
        2,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "({}) Observed: {}   expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e2r.RabinowitzRichter[1]()])
def test_show(scheme):
    quadpy.e2r.show(scheme)
    return


if __name__ == "__main__":
    scheme_ = quadpy.e2r.RabinowitzRichter(5)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
