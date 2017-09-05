# -*- coding: utf-8 -*-
#
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.enr.Stroud(n, key), 1.0e-14) for n in range(4, 6)
     for key in quadpy.enr.Stroud.keys
     ]
    + [(quadpy.enr.StroudSecrest(n, key), 1.0e-14) for n in range(2, 6)
       for key in quadpy.enr.StroudSecrest.keys
       ]
    )
def test_scheme(scheme, tol):
    n = scheme.dim
    degree = check_degree(
            lambda poly: quadpy.enr.integrate(poly, scheme),
            integrate_monomial_over_enr,
            lambda k: quadpy.helpers.partition(k, n),
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    dim_ = 2
    # quadpy.e3r2.show(quadpy.enr.Stroud(dim_, '5-1a'), backend='vtk')
    scheme_ = quadpy.enr.Stroud(dim_, '5-3')
    test_scheme(scheme_, 1.0e-14)
