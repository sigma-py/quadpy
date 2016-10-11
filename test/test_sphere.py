# -*- coding: utf-8 -*-
#
import math
import numpy
import numpy.testing
import pytest
import quadrature
import sympy

from test_tetrahedron import _create_monomial_exponents

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integral_monomial_over_unit_sphere(alpha):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448.
    '''
    if any(alpha % 2 == 1):
        return 0.0

    def gamma2(kk):
        '''Gamma function with the half argument, i.e.,
        gamma2(k) = gamma(k/2)
        '''
        if kk % 2 == 0:
            return math.factorial(kk/2 - 1)

        k = int(kk / 2.0)
        return numpy.sqrt(numpy.pi) \
            * numpy.prod([(i + 0.5) for i in range(k)])

    num = 2 * numpy.prod([gamma2(a + 1) for a in alpha])

    return num / gamma2(sum(alpha) + len(alpha))


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


@pytest.mark.parametrize('scheme', [
    quadrature.sphere.Lebedev(1),
    quadrature.sphere.Lebedev(2),
    quadrature.sphere.Lebedev(3),
    quadrature.sphere.Lebedev(4),
    ])
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    midpoint = numpy.array([0.0, 0.0, 0.0])
    radius = 1.0
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in _create_monomial_exponents(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1] * x[2]**k[2]
            # exact_val = _integrate_exact(poly, midpoint, radius)
            exact_val = _integral_monomial_over_unit_sphere(k)
            val = quadrature.sphere.integrate(
                    poly, midpoint, radius, scheme
                    )
            if abs(exact_val - val) > 1.0e-10:
                success = False
                break
        if not success:
            break
        if degree >= max_degree:
            break
        degree += 1
    assert degree >= max_degree
    return


def test_show():
    quadrature.sphere.show(
        quadrature.sphere.Lebedev(4)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
