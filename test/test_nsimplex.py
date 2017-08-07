# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy

from helpers import (
    check_degree, integrate_monomial_over_standard_simplex
    )


@pytest.mark.parametrize(
    'scheme',
    [quadpy.nsimplex.GrundmannMoeller(3, k) for k in range(5)] +
    [quadpy.nsimplex.GrundmannMoeller(4, k) for k in range(5)] +
    [quadpy.nsimplex.GrundmannMoeller(5, k) for k in range(5)] +
    [quadpy.nsimplex.GrundmannMoeller(6, k) for k in range(5)] +
    #
    [quadpy.nsimplex.Walkington(3, k) for k in [1, 2, 3, 5, 7]] +
    [quadpy.nsimplex.Walkington(4, k) for k in [1, 2, 3]] +
    [quadpy.nsimplex.Walkington(5, k) for k in [1, 2, 3]] +
    [quadpy.nsimplex.Walkington(6, k) for k in [1, 2, 3]]
    )
def test_scheme(scheme):
    n = scheme.dim
    simplex = numpy.zeros((n+1, n))
    for k in range(n):
        simplex[k+1, k] = 1.0
    degree = check_degree(
            lambda poly: quadpy.nsimplex.integrate(poly, simplex, scheme),
            integrate_monomial_over_standard_simplex,
            lambda k: quadpy.helpers.partition(k, n),
            scheme.degree + 1
            )
    assert degree >= scheme.degree
    return


if __name__ == '__main__':
    n_ = 3
    scheme_ = quadpy.nsimplex.GrundmannMoeller(n_, 5)
    test_scheme(scheme_)
