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


@pytest.mark.parametrize('scheme', [
    quadrature.line_segment.Midpoint(),
    quadrature.line_segment.Trapezoidal(),
    quadrature.line_segment.GaussLegendre(1),
    quadrature.line_segment.GaussLegendre(2),
    quadrature.line_segment.GaussLegendre(3),
    quadrature.line_segment.GaussLegendre(4),
    quadrature.line_segment.GaussLegendre(5),
    quadrature.line_segment.GaussPatterson(0),
    quadrature.line_segment.GaussPatterson(1),
    quadrature.line_segment.GaussPatterson(2),
    quadrature.line_segment.GaussPatterson(3),
    quadrature.line_segment.GaussPatterson(4),
    quadrature.line_segment.GaussPatterson(5),
    quadrature.line_segment.GaussPatterson(6),
    quadrature.line_segment.ClenshawCurtis(1),
    quadrature.line_segment.ClenshawCurtis(2),
    quadrature.line_segment.ClenshawCurtis(3),
    quadrature.line_segment.ClenshawCurtis(4),
    quadrature.line_segment.ClenshawCurtis(5),
    quadrature.line_segment.ClenshawCurtis(9),
    quadrature.line_segment.ClenshawCurtis(17),
    quadrature.line_segment.ClenshawCurtis(33),
    quadrature.line_segment.ClenshawCurtis(65),
    quadrature.line_segment.NewtonCotesClosed(1),
    quadrature.line_segment.NewtonCotesClosed(2),
    quadrature.line_segment.NewtonCotesClosed(3),
    quadrature.line_segment.NewtonCotesClosed(4),
    quadrature.line_segment.NewtonCotesOpen(1),
    quadrature.line_segment.NewtonCotesOpen(2),
    quadrature.line_segment.NewtonCotesOpen(3),
    quadrature.line_segment.NewtonCotesOpen(4),
    ])
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


def test_show():
    quadrature.line_segment.show(
        0.0, 1.0,
        # quadrature.line_segment.NewtonCotesOpen(6),
        # quadrature.line_segment.NewtonCotesClosed(15),
        # quadrature.line_segment.GaussLegendre(31),
        # quadrature.line_segment.GaussPatterson(4),
        quadrature.line_segment.ClenshawCurtis(33),
        render=False
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
    # scheme = quadrature.line_segment.GaussLegendre(7)
    # check_scheme(scheme)
