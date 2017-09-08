# -*- coding: utf-8 -*-
#
import math

import pytest
import quadpy

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.e1r.GaussLaguerre(n), 1.0e-14) for n in range(1, 10)]
    )
def test_scheme(scheme, tol):
    degree = check_degree(
            lambda poly: quadpy.e1r.integrate(poly, scheme),
            lambda k: math.factorial(k[0]),
            lambda k: quadpy.helpers.partition(k, 1),
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


if __name__ == '__main__':
    scheme_ = quadpy.e1r.GaussLaguerre(5)
    test_scheme(scheme_, 1.0e-14)
