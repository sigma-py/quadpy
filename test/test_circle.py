# -*- coding: utf-8 -*-
#
import pytest
import quadpy

from helpers import integrate_monomial_over_unit_sphere, check_degree


@pytest.mark.parametrize(
    'scheme',
    [quadpy.circle.Equidistant(k) for k in range(1, 6)]
    )
def test_scheme(scheme):
    degree = check_degree(
            lambda poly: quadpy.circle.integrate(
                poly, [0.0, 0.0], 1.0, scheme
                ),
            integrate_monomial_over_unit_sphere,
            lambda n: quadpy.helpers.partition(n, 2),
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.circle.Equidistant(3)]
    )
def test_show(scheme):
    quadpy.circle.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.circle.Equidistant(30)
    test_scheme(scheme_)
    test_show(scheme_)
