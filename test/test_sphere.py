# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3, check_degree
import math
from matplotlib import pyplot as plt
import numpy
import pytest
import quadpy
import sympy


def _integral_monomial_over_unit_sphere(alpha):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448.
    '''
    if any(alpha % 2 == 1):
        return 0.0

    # Use lgamma since other with ordinary gamma, numerator and denominator
    # might overflow.
    return 2.0 * math.exp(
        math.fsum([math.lgamma(0.5*(a+1)) for a in alpha])
        - math.lgamma(math.fsum([0.5*(a+1) for a in alpha]))
        )


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
    [quadpy.sphere.Lebedev(k) for k in range(1, 33)]
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
            _integral_monomial_over_unit_sphere,
            create_monomial_exponents3,
            min(30, scheme.degree + 1)
            )
    assert degree >= min(30, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(4)]
    )
def test_show(scheme):
    quadpy.sphere.show(scheme)
    return


if __name__ == '__main__':
    scheme = quadpy.sphere.Lebedev(4)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
