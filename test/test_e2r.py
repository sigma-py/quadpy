# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.helpers import kahan_dot

from helpers import check_degree, integrate_monomial_over_enr


@pytest.mark.parametrize(
    'scheme,tol',
    [(quadpy.e2r.RabinowitzRichter(k), 1.0e-14)
     for k in quadpy.e2r.RabinowitzRichter.keys
     ]
    + [(quadpy.e2r.Stroud(key), 1.0e-14) for key in quadpy.e2r.Stroud.keys]
    + [(quadpy.e2r.StroudSecrest(key), 1.0e-14)
       for key in quadpy.e2r.StroudSecrest.keys
       ]
    )
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
            lambda poly: quadpy.e2r.integrate(poly, scheme, dot=kahan_dot),
            integrate_monomial_over_enr,
            2,
            scheme.degree + 1,
            tol=tol
            )
    assert degree == scheme.degree, \
        'Observed: {}   expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.e2r.RabinowitzRichter(1)]
    )
def test_show(scheme):
    quadpy.e2r.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.e2r.RabinowitzRichter(5)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)
