# -*- coding: utf-8 -*-
#
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
        sympy.integrate(rtheta_x_rphi * f(x_xi), (phi, 0.0, numpy.pi)),
        (theta, 0, 2*numpy.pi)
        )
    return float(exact)


@pytest.mark.parametrize('scheme', [
    quadrature.sphere.Lebedev(1),
    quadrature.sphere.Lebedev(2),
    quadrature.sphere.Lebedev(3),
    quadrature.sphere.Lebedev(4),
    quadrature.sphere.Lebedev(5),
    quadrature.sphere.Lebedev(6),
    quadrature.sphere.Lebedev(7),
    quadrature.sphere.Lebedev(8),
    quadrature.sphere.Lebedev(9),
    quadrature.sphere.Lebedev(10),
    quadrature.sphere.Lebedev(11),
    quadrature.sphere.Lebedev(12),
    quadrature.sphere.Lebedev(13),
    quadrature.sphere.Lebedev(14),
    quadrature.sphere.Lebedev(15),
    quadrature.sphere.Lebedev(16),
    quadrature.sphere.Lebedev(17),
    quadrature.sphere.Lebedev(18),
    quadrature.sphere.Lebedev(19),
    quadrature.sphere.Lebedev(20),
    quadrature.sphere.Lebedev(21),
    quadrature.sphere.Lebedev(22),
    quadrature.sphere.Lebedev(23),
    quadrature.sphere.Lebedev(24),
    quadrature.sphere.Lebedev(25),
    quadrature.sphere.Lebedev(26),
    quadrature.sphere.Lebedev(27),
    quadrature.sphere.Lebedev(28),
    quadrature.sphere.Lebedev(29),
    quadrature.sphere.Lebedev(30),
    quadrature.sphere.Lebedev(31),
    quadrature.sphere.Lebedev(32),
    ])
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    midpoint = numpy.array([0.0, 0.0, 0.0])
    radius = 1.0
    success = True
    degree = 0
    # Cap the tested degree. That's ugly, but otherwise tests just take
    # forever.
    # TODO come up with a better solution here
    max_degree = min(10, scheme.degree + 1)
    while success:
        for k in _create_monomial_exponents(degree):
            def poly(x):
                return x[0]**k[0] + x[1]**k[1] + x[2]**k[2]
            print(k)
            exact_val = _integrate_exact(poly, midpoint, radius)
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
    # numpy.testing.assert_equal(degree-1, scheme.degree)
    return


def test_show():
    quadrature.sphere.show(
        quadrature.sphere.Lebedev(13)
        )
    return


if __name__ == '__main__':
    # test_show()
    # plt.show()
    scheme = quadrature.sphere.Lebedev(1)
    test_scheme(scheme)
