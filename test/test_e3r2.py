# -*- coding: utf-8 -*-
#
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr2


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.e3r2.Stroud(index), 1.0e-14) for index in [
        '5-1', '5-2a', '5-2b', '5-3',
        '7-1a', '7-1b', '7-2a', '7-2b',
        '14-1',
        ]]
    )
def test_scheme(scheme, tol):
    degree = check_degree(
            lambda poly: quadpy.e3r2.integrate(poly, scheme),
            integrate_monomial_over_enr2,
            lambda n: quadpy.helpers.partition(n, 3),
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.e3r2.Stroud('5-1')]
    )
def test_show(scheme, backend='mpl'):
    quadpy.e3r2.show(scheme, backend=backend)
    return


if __name__ == '__main__':
    scheme_ = quadpy.e3r2.Stroud('7-2b')
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_, backend='vtk')
