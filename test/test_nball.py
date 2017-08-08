# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nball

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme',
    [quadpy.nball.Dobrodeev1970(n) for n in range(3, 10)]
    + [quadpy.nball.Dobrodeev1978(n) for n in range(2, 7)]
    + [quadpy.nball.Stroud(dim, index) for dim in range(2, 10) for index in [
        'Sn 2-1',
        'Sn 3-1', 'Sn 3-2',
        'Sn 5-2', 'Sn 5-3', 'Sn 5-4', 'Sn 5-5', 'Sn 5-6',
        ]]
    + [quadpy.nball.Stroud(dim, index) for dim in range(4, 8) for index in [
        'Sn 5-1a', 'Sn 5-1b',
        ]]
    )
def test_scheme(scheme):
    tol = 1.0e-14
    n = scheme.dim
    center = numpy.zeros(n)
    radius = 1.0
    degree = check_degree(
            lambda poly: quadpy.nball.integrate(poly, center, radius, scheme),
            integrate_monomial_over_unit_nball,
            lambda k: quadpy.helpers.partition(k, n),
            scheme.degree + 1,
            tol=tol
            )
    assert degree >= scheme.degree, \
        'observed: {}, expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    n_ = 8
    scheme_ = quadpy.nball.Stroud(n_, 'Sn 5-5')
    test_scheme(scheme_)
