# -*- coding: utf-8 -*-
#
from helpers import \
        create_monomial_exponents2, \
        integrate_monomial_over_unit_circle

import numpy
import numpy.testing
import pytest
import quadrature

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integrate_exact(k):
    '''We have

    I = \int_0^1 \int_0^2pi r * (r cos(phi))**k[0] (r sin(phi))**k[1]
      = 1.0/(2+k[0]+k[1]) * \int_0^2pi cos(phi)**k[0] sin(phi)**k[1]
    '''
    return 1.0/(2 + k[0] + k[1]) * integrate_monomial_over_unit_circle(k)


@pytest.mark.parametrize(
    'scheme',
    [quadrature.disk.Peirce(k) for k in range(1, 6)]
    + [quadrature.disk.Lether(k) for k in range(1, 6)]
    )
def test_scheme(scheme):
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in create_monomial_exponents2(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1]
            exact_val = _integrate_exact(k)
            val = quadrature.disk.integrate(poly, scheme)
            # print('k, exact_val', k, exact_val, val)
            if abs(exact_val - val) > 1.0e-10:
                success = False
                break
        if not success:
            break
        if degree >= max_degree:
            break
        degree += 1
    numpy.testing.assert_equal(degree-1, scheme.degree)
    return


def test_show():
    quadrature.disk.show(
        quadrature.disk.Peirce(3)
        # quadrature.disk.Lether(5)
        )
    return

if __name__ == '__main__':
    test_show()
    plt.show()
    # scheme = From1d(quadrature.line_segment.NewtonCotesClosed(15))
    scheme = quadrature.disk.Lether(5)
    test_scheme(scheme)
