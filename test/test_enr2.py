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
        '5-1a',
        '5-2',
        ]]
    + [(quadpy.enr2.Stroud(n, index), 1.0e-14) for n in [3, 5, 6]
        for index in [
        '5-1b',
        ]]
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
    dim_ = 7
    # quadpy.e3r2.show(quadpy.enr2.Stroud(dim_, '5-1a'), backend='vtk')
    scheme_ = quadpy.enr2.Stroud(dim_, '5-1a')
    test_scheme(scheme_, 1.0e-14)
