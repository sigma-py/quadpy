# -*- coding: utf-8 -*-
#
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nball

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.ball.HammerStroud(k), 1.0e-14) for k in [
        '11-3', '12-3', '14-3a', '14-3b', '15-3a', '15-3b',
        ]]
    + [(quadpy.ball.Stroud(k), 1.0e-14) for k in [
        'S3 3-1',
        'S3 5-1', 'S3 5-2',
        'S3 7-1a', 'S3 7-1b', 'S3 7-2', 'S3 7-3', 'S3 7-4',
        'S3 14-1',
        ]]
    )
def test_scheme(scheme, tol):
    degree = check_degree(
            lambda poly: quadpy.ball.integrate(
                poly, [0.0, 0.0, 0.0], 1.0, scheme
                ),
            integrate_monomial_over_unit_nball,
            3,
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
def test_show(scheme, backend='mpl'):
    quadpy.ball.show(scheme, backend=backend)
    return


if __name__ == '__main__':
    scheme_ = quadpy.ball.Stroud('S3 14-1', symbolic=True)
    test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_, backend='vtk')
