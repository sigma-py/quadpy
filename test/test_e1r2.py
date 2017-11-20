# -*- coding: utf-8 -*-
#
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr2


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.e1r2.GaussHermite(n), 1.0e-14) for n in range(1, 10)]
    )
def test_scheme(scheme, tol):
    degree = check_degree(
            lambda poly: quadpy.e1r2.integrate(poly, scheme),
            integrate_monomial_over_enr2,
            lambda k: quadpy.helpers.partition(k, 1),
            1,
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.e1r2.GaussHermite(2)]
    )
def test_show(scheme):
    quadpy.e1r2.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.e1r2.GaussHermite(10)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
