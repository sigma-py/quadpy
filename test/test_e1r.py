# -*- coding: utf-8 -*-
#
import math

import numpy
import pytest
import quadpy

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.e1r.GaussLaguerre(n), 1.0e-14) for n in range(1, 10)]
    )
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
            lambda poly: quadpy.e1r.integrate(poly, scheme),
            lambda k: math.factorial(k[0]),
            1,
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.e1r.GaussLaguerre(1)]
    )
def test_show(scheme):
    quadpy.e1r.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.e1r.GaussLaguerre(3)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
