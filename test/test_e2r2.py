# -*- coding: utf-8 -*-
#
import math
import numpy
import pytest
import quadpy

from helpers import check_degree


def integrate_monomial_over_enr2(k):
    if numpy.any(k % 2 == 1):
        return 0
    return numpy.prod([math.gamma((kk+1) / 2.0) for kk in k])


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.e2r2.RabinowitzRichter(k), 1.0e-14) for k in range(1, 2)]
    )
def test_scheme(scheme, tol):
    degree = check_degree(
            lambda poly: quadpy.e2r2.integrate(poly, scheme),
            integrate_monomial_over_enr2,
            lambda n: quadpy.helpers.partition(n, 2),
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.e2r2.RabinowitzRichter(1)]
    )
def test_show(scheme):
    quadpy.e2r2.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.e2r2.RabinowitzRichter(1)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
