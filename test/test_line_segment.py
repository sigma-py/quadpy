# -*- coding: utf-8 -*-
#
import quadrature

import os
import matplotlib as mpl
import pytest
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.Midpoint()]
    + [quadrature.line_segment.Trapezoidal()]
    + [quadrature.line_segment.GaussLegendre(k) for k in range(1, 6)]
    + [quadrature.line_segment.GaussPatterson(0) for k in range(7)]
    + [quadrature.line_segment.ClenshawCurtis(k) for k in [
        1, 2, 3, 4, 5, 9, 17, 33, 65
        ]]
    + [quadrature.line_segment.NewtonCotesClosed(k) for k in range(1, 5)]
    + [quadrature.line_segment.NewtonCotesOpen(k) for k in range(1, 5)]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    degree = 0
    while True:
        def poly(x): return x**degree
        # Set bounds such that the values are between 0.5 and 1.5.
        a = 0.5**(1.0/(degree+1))
        b = 1.5**(1.0/(degree+1))
        exact_val = 1.0/(degree+1)
        val = quadrature.line_segment.integrate(poly, a, b, scheme)
        if abs(exact_val - val) / abs(exact_val) > 1.0e-12:
            break
        if degree >= scheme.degree:
            break
        degree += 1
    assert degree >= scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadrature.line_segment.NewtonCotesClosed(5)]
    )
def test_show(scheme):
    quadrature.line_segment.show(0.0, 1.0, scheme)
    return


if __name__ == '__main__':
    scheme = quadrature.line_segment.GaussLegendre(7)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
