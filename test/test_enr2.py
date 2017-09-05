# -*- coding: utf-8 -*-
#
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr2


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.enr2.Stroud(n, index), 1.0e-14) for n in range(2, 8)
        for index in [
        '3-1', '3-2',
        '5-1a', '5-2', '5-4',
        '5-5a',
        ]]
    + [(quadpy.enr2.Stroud(n, '5-3'), 1.0e-14) for n in range(3, 8)]
    + [(quadpy.enr2.Stroud(n, '5-5b'), 1.0e-14) for n in [2]]
    + [(quadpy.enr2.Stroud(n, '5-6'), 1.0e-14) for n in range(5, 8)]
    + [(quadpy.enr2.Stroud(n, '5-1b'), 1.0e-14) for n in [3, 5, 6]]
    + [(quadpy.enr2.Stroud(n, '7-1a'), 1.0e-14) for n in [2, 3, 4, 6, 7]]
    + [(quadpy.enr2.Stroud(n, '7-1b'), 1.0e-14) for n in [3, 4]]
    + [(quadpy.enr2.Stroud(n, '7-2'), 1.0e-14) for n in [3, 4, 5, 6]]
    + [(quadpy.enr2.Stroud(n, index), 1.0e-14) for n in [3, 4, 5, 6]
        for index in ['7-3b', '9-1a']]
    + [(quadpy.enr2.Stroud(n, '7-3a'), 1.0e-14) for n in [3, 4, 6]]
    + [(quadpy.enr2.Stroud(n, '9-1b'), 1.0e-14) for n in [4, 5, 6]]
    + [(quadpy.enr2.Stroud(n, '11-1a'), 1.0e-14) for n in [3, 4]]
    + [(quadpy.enr2.Stroud(n, '11-1b'), 1.0e-14) for n in [3, 4, 5]]
    + [(quadpy.enr2.StroudSecrest(n, 'II'), 1.0e-14) for n in range(2, 6)]
    )
def test_scheme(scheme, tol):
    n = scheme.dim
    degree = check_degree(
            lambda poly: quadpy.enr2.integrate(poly, scheme),
            integrate_monomial_over_enr2,
            lambda k: quadpy.helpers.partition(k, n),
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    dim_ = 5
    # quadpy.e3r2.show(quadpy.enr2.Stroud(dim_, '5-1a'), backend='vtk')
    scheme_ = quadpy.enr2.Stroud(dim_, '11-1b')
    test_scheme(scheme_, 1.0e-14)
