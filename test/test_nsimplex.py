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
    [quadpy.nsimplex.GrundmannMoeller(dim, k)
     for dim in range(3, 7)
     for k in range(5)
     ]
    #
    + [quadpy.nsimplex.Stroud(dim, index)
       for dim in range(3, 7)
       for index in ['Tn 1-1']
       ]
    #
    + [quadpy.nsimplex.Walkington(3, k) for k in [1, 2, 3, 5, 7]]
    + [quadpy.nsimplex.Walkington(4, k)
       for dim in range(4, 7)
       for k in [1, 2, 3]
       ]
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
    assert degree >= scheme.degree, \
        'Observed: {}, expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    n_ = 4
    scheme_ = quadpy.nsimplex.Stroud(n_, 'Tn 3-4')
    test_scheme(scheme_)
