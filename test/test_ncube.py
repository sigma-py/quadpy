# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy

from helpers import check_degree


stroud_idx = [
    'Cn 1-1', 'Cn 1-2',
    'Cn 2-1', 'Cn 2-2',
    'Cn 3-1',
    'Cn 3-2', 'Cn 3-3', 'Cn 3-4', 'Cn 3-5', 'Cn 3-6',
    'Cn 5-2', 'Cn 5-3', 'Cn 5-4', 'Cn 5-5', 'Cn 5-6', 'Cn 5-7', 'Cn 5-8',
    'Cn 5-9',
    # 'Cn 7-1',
    ]


@pytest.mark.parametrize(
    'scheme',
    [quadpy.ncube.Stroud(3, k) for k in stroud_idx]
    + [quadpy.ncube.Stroud(4, k) for k in stroud_idx]
    + [quadpy.ncube.Stroud(5, k) for k in stroud_idx]
    + [quadpy.ncube.Stroud(6, k) for k in stroud_idx]
    )
def test_scheme(scheme):
    n = scheme.dim
    ncube_limits = [[0.0, 1.0]] * n
    ncube = quadpy.ncube.ncube_points(*ncube_limits)
    degree = check_degree(
            lambda poly: quadpy.ncube.integrate(poly, ncube, scheme),
            lambda exp: _integrate_monomial_over_ncube(ncube_limits, exp),
            lambda k: quadpy.helpers.partition(k, n),
            scheme.degree + 1
            )
    assert degree >= scheme.degree, \
        'observed: {}, expected: {}'.format(degree, scheme.degree)
    return


def _integrate_monomial_over_ncube(ncube_limits, exp):
    return numpy.prod([
            (a[1]**(k+1) - a[0]**(k+1)) / (k+1) for a, k in
            zip(ncube_limits, exp)
            ])


if __name__ == '__main__':
    n = 5
    scheme_ = quadpy.ncube.Stroud(n, 'Cn 3-4')
    test_scheme(scheme_)
