# -*- coding: utf-8 -*-
#
from helpers import \
        create_monomial_exponents2, \
        integrate_monomial_over_unit_circle, \
        check_degree

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
    degree = check_degree(
            lambda poly: quadrature.circle.integrate(poly, scheme),
            integrate_monomial_over_unit_circle,
            scheme.degree + 1
            )
    numpy.testing.assert_equal(degree, scheme.degree)
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
