# -*- coding: utf-8 -*-
#
import pytest
import quadpy
from quadpy.ncube.helpers import integrate_monomial_over_ncube

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme, tol',
    [(quadpy.ncube.Dobrodeev1970(dim), 1.0e-14) for dim in range(5, 8)]
    + [(quadpy.ncube.Dobrodeev1978(dim), 1.0e-14) for dim in range(2, 8)]
    + [(quadpy.ncube.HammerStroud(n, k), 1.0e-14) for n in range(3, 7)
       for k in ['1-n', '2-n']
       ]
    + [(quadpy.ncube.Stroud(n, k), 1.0e-14) for n in range(3, 7)
       for k in [
        'Cn 1-1', 'Cn 1-2',
        'Cn 2-1', 'Cn 2-2',
        'Cn 3-1',
        'Cn 3-2', 'Cn 3-3', 'Cn 3-4', 'Cn 3-5', 'Cn 3-6',
        'Cn 5-2', 'Cn 5-3', 'Cn 5-4', 'Cn 5-5', 'Cn 5-6', 'Cn 5-7', 'Cn 5-8',
        'Cn 5-9',
        'Cn 7-1',
        ]
       ]
    )
def test_scheme(scheme, tol):
    n = scheme.dim
    ncube_limits = [[0.0, 1.0]] * n
    ncube = quadpy.ncube.ncube_points(*ncube_limits)
    degree = check_degree(
            lambda poly: quadpy.ncube.integrate(poly, ncube, scheme),
            lambda exp: integrate_monomial_over_ncube(ncube_limits, exp),
            lambda k: quadpy.helpers.partition(k, n),
            n,
            scheme.degree + 1,
            tol=tol
            )
    assert degree >= scheme.degree, \
        'observed: {}, expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    n_ = 4
    scheme_ = quadpy.ncube.Stroud(n_, 'Cn 7-1')
    test_scheme(scheme_, 1.0e-14)
