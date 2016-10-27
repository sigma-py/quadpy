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


@pytest.mark.parametrize(
    'scheme',
    [quadrature.circle.Equidistant(k) for k in range(1, 6)]
    )
def test_scheme(scheme):
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in create_monomial_exponents2(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1]
            exact_val = integrate_monomial_over_unit_circle(k)
            val = quadrature.circle.integrate(poly, scheme)
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


@pytest.mark.parametrize(
    'scheme',
    [quadrature.circle.Equidistant(3)]
    )
def test_show(scheme):
    quadrature.circle.show(scheme)
    return

if __name__ == '__main__':
    scheme = quadrature.circle.Equidistant(30)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
