# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nsphere

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme',
    [quadpy.circle.Krylov(k, symbolic=False) for k in range(1, 6)]
    )
def test_scheme(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
            lambda poly: quadpy.circle.integrate(
                poly, [0.0, 0.0], 1.0, scheme
                ),
            integrate_monomial_over_unit_nsphere,
            2,
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.circle.Krylov(3)]
    )
def test_show(scheme):
    quadpy.circle.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.circle.Krylov(30)
    test_scheme(scheme_)
    test_show(scheme_)
