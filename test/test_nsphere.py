# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.nsphere.helpers import integrate_monomial_over_unit_nsphere

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme',
    [quadpy.nsphere.Dobrodeev1978(n) for n in range(2, 7)]
    + [quadpy.nsphere.Stroud(n, index) for n in range(2, 7) for index in [
        'Un 3-1',
        ]]
    + [quadpy.nsphere.Stroud1967(n) for n in range(2, 7)]
    )
def test_scheme(scheme):
    tol = 1.0e-14
    n = scheme.dim
    center = numpy.zeros(n)
    rad = 1.0
    degree = check_degree(
            lambda poly: quadpy.nsphere.integrate(poly, center, rad, scheme),
            integrate_monomial_over_unit_nsphere,
            lambda k: quadpy.helpers.partition(k, n),
            scheme.degree + 1,
            tol=tol
            )
    assert degree >= scheme.degree, \
        'observed: {}, expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    n_ = 5
    scheme_ = quadpy.nsphere.Stroud(n_, 'Un 3-1')
    test_scheme(scheme_)
