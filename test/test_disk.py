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
    degree = check_degree(
            lambda poly: quadrature.disk.integrate(poly, scheme),
            _integrate_exact,
            scheme.degree + 1
            )
    numpy.testing.assert_equal(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.disk.Lether(3)]
    )
def test_show(scheme):
    quadrature.disk.show(scheme)
    return

if __name__ == '__main__':
    scheme = quadrature.disk.Lether(5)
    # scheme = From1d(quadrature.line_segment.NewtonCotesClosed(15))
    test_scheme(scheme)
    # test_show(scheme)
    # plt.show()
