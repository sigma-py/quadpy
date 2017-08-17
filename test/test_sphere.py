# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nsphere
import sympy

from helpers import check_degree


def _integrate_exact(f, midpoint, radius):
    #
    # New coordinates: phi, theta.
    #
    phi = sympy.Symbol('phi')
    theta = sympy.Symbol('theta')
    x_xi = [
        midpoint[0] + radius*sympy.sin(phi)*sympy.cos(theta),
        midpoint[1] + radius*sympy.sin(phi)*sympy.sin(theta),
        midpoint[2] + radius*sympy.cos(phi),
        ]
    rtheta_x_rphi = sympy.sin(phi) * radius**2
    exact = sympy.integrate(
        sympy.integrate(rtheta_x_rphi * f(x_xi), (phi, 0.0, sympy.pi)),
        (theta, 0, 2*sympy.pi)
        )
    return float(exact)


@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(degree) for degree in [
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,
        59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131
        ]]
    + [quadpy.sphere.McLaren()]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    midpoint = numpy.array([0.0, 0.0, 0.0])
    radius = 1.0
    degree = check_degree(
            lambda poly: quadpy.sphere.integrate(
                poly, midpoint, radius, scheme, sumfun=numpy.sum
                ),
            integrate_monomial_over_unit_nsphere,
            lambda n: quadpy.helpers.partition(n, 3),
            min(30, scheme.degree + 1)
            )
    assert degree >= min(30, scheme.degree), \
        'Observed: {}, expected: {}'.format(degree, scheme.degree)
    return


# some basic sanity tests for integration_spherical
@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(degree) for degree in [
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,
        59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131
        ]]
    )
def test_scheme_spherical(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    res = quadpy.sphere.integrate_spherical(
            lambda phi_theta: 1.0,
            radius=1.0,
            rule=scheme
            )
    assert abs(res - 4*numpy.pi) < 1.0e-11
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(7)]
    )
def test_show(scheme):
    quadpy.sphere.show(scheme)
    return


if __name__ == '__main__':
    scheme_ = quadpy.sphere.McLaren()
    test_scheme(scheme_)
    # test_scheme_spherical(scheme_)
    test_show(scheme_)
