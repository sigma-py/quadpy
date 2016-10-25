# -*- coding: utf-8 -*-
#
import math
import numpy
import numpy.testing
import pytest
import quadrature
import sympy

from test_triangle import _create_monomial_exponents

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integrate_exact(k):
    '''The integral

    I = \int_0^1 \int_0^2pi r * (r cos(phi))**k[0] (r sin(phi))**k[1]

    equals 0 if any of k[0], k[1] is odd. If both are even, then

    I = 1.0/(2+k[0]+k[1]) * \int_0^2pi sin**2m * cos**2n
      = 1.0/(2+k[0]+k[1]) * \int_0^2pi (1-cos**2)**m * cos**2n

    The latter integral is
        \int_0^2pi (1-cos**2)**m * cos**2n
      = \sum_k=0^m (m over k) \int cos**2(k+n).

    With
      \int cos**n = (n-1)/n \int cos**(n-2)
    one gets
        \int_0^2pi cos**n
      = (n-1)!! / n!! * 2*pi
      = (2n)! / (2**n * n!)**2 * 2*pi,
    so
      I = 1.0/(2 + k[0] + k[1]) \
        * \sum_k=0^m (-1)^i(m over k) (2(k+n))! / (2**(k+n) * (k+n)!)**2 * 2*pi

    The quotient can be computed via log-gamma, i.e.,

      (2d)! / (2**d * d!)**2
    = exp(lgamma(2d) - 2d log(2) - 2 lgamma(d)).
    '''
    if any(numpy.array(k) % 2 == 1):
        return 0.0

    m = k[0] // 2
    n = k[1] // 2
    return 1.0/(2 + k[0] + k[1]) \
        * 2*math.pi * sum([
            (-1)**i
            * sympy.binomial(m, i)
            * math.factorial(2*(i+n))
            / (2**(i+n) * math.factorial(i+n))**2
            for i in range(m+1)
            ])
    # return 1.0/(2 + k[0] + k[1]) \
    #     * 2*math.pi * sum([
    #         (-1)**i
    #         * sympy.binomial(m, i)
    #         * math.exp(
    #             + math.lgamma(2*(i+n))
    #             - 2*math.lgamma(i+n)
    #             - 2*(i+n)*math.log(2.0)
    #             )
    #         ] for i in range(m+1)
    #         )


@pytest.mark.parametrize('scheme', [
    quadrature.circle.Peirce(1),
    quadrature.circle.Peirce(2),
    quadrature.circle.Peirce(3),
    quadrature.circle.Peirce(4),
    quadrature.circle.Peirce(5),
    quadrature.circle.Lether(1),
    quadrature.circle.Lether(2),
    quadrature.circle.Lether(3),
    quadrature.circle.Lether(4),
    quadrature.circle.Lether(5),
    ])
def test_scheme(scheme):
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in _create_monomial_exponents(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1]
            exact_val = _integrate_exact(k)
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


def test_show():
    quadrature.circle.show(
        quadrature.circle.Lether(5)
        )
    return

if __name__ == '__main__':
    test_show()
    plt.show()
    # scheme = From1d(quadrature.line.NewtonCotesClosed(15))
    scheme = quadrature.circle.Lether(5)
    test_scheme(scheme)
