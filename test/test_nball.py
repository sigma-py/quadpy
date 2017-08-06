# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_unit_nball


@pytest.mark.parametrize(
    'scheme',
    [quadpy.nball.Dobrodeev(n) for n in range(3, 10)]
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
    n_ = 3
    scheme_ = quadpy.nball.Dobrodeev(n_)
    test_scheme(scheme_)
