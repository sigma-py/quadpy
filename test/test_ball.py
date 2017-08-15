# -*- coding: utf-8 -*-
#
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nball

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.ball.HammerStroud(k), 1.0e-14) for k in [
        '11-3', '12-3', '14-3a', '14-3b',
        ]]
    )
def test_scheme(scheme, tol):
    degree = check_degree(
            lambda poly: quadpy.ball.integrate(
                poly, [0.0, 0.0, 0.0], 1.0, scheme
                ),
            integrate_monomial_over_unit_nball,
            lambda n: quadpy.helpers.partition(n, 3),
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.ball.HammerStroud('11-3')]
    )
def test_show(scheme):
    quadpy.ball.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.ball.HammerStroud('14-3b')
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
