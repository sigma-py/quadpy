# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3
import math
import numpy
import numpy.testing
import pytest
import quadrature
import sympy

import os
import matplotlib as mpl
if 'DISPLAY' not in os.environ:
    # headless mode, for remote executions (and travis)
    mpl.use('Agg')
from matplotlib import pyplot as plt


def _integrate_exact(f, tetrahedron):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference tetrahedron [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the tetrahedron. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    xi = sympy.DeferredVector('xi')
    x_xi = \
        + tetrahedron[0] * (1 - xi[0] - xi[1] - xi[2]) \
        + tetrahedron[1] * xi[0] \
        + tetrahedron[2] * xi[1] \
        + tetrahedron[3] * xi[2]
    abs_det_J = 6 * quadrature.tetrahedron.volume(tetrahedron)
    exact = sympy.integrate(
        sympy.integrate(
          sympy.integrate(abs_det_J * f(x_xi), (xi[2], 0, 1-xi[0]-xi[1])),
          (xi[1], 0, 1-xi[0])
          ),
        (xi[0], 0, 1)
      )
    return float(exact)


def _integrate_monomial_over_standard_tet(k):
    '''The integral of monomials over the standard tetrahedron is given by

    \int_T x_0^k0 * x1^k1 * x2^k2 = (k0!*k1!*k2!) / (3+k0+k1+k2)!,

    see, e.g.,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Linbo Zhang, Tao Cui and Hui Liu,
    Journal of Computational Mathematics,
    Vol. 27, No. 1 (January 2009), pp. 89-96
    '''
    # exp-log to account for large values in numerator and denominator
    return math.exp(
        math.fsum([math.lgamma(kk+1) for kk in k])
        - math.lgamma(4 + sum(k))
        )


@pytest.mark.parametrize('scheme', [
    quadrature.tetrahedron.Keast(0),
    quadrature.tetrahedron.Keast(1),
    quadrature.tetrahedron.Keast(2),
    quadrature.tetrahedron.Keast(3),
    quadrature.tetrahedron.Keast(4),
    quadrature.tetrahedron.Keast(5),
    quadrature.tetrahedron.Keast(6),
    quadrature.tetrahedron.Keast(7),
    quadrature.tetrahedron.Keast(8),
    quadrature.tetrahedron.Keast(9),
    quadrature.tetrahedron.Keast(10),
    quadrature.tetrahedron.NewtonCotesClosed(1),
    quadrature.tetrahedron.NewtonCotesClosed(2),
    quadrature.tetrahedron.NewtonCotesClosed(3),
    quadrature.tetrahedron.NewtonCotesClosed(4),
    quadrature.tetrahedron.NewtonCotesClosed(5),
    quadrature.tetrahedron.NewtonCotesClosed(6),
    quadrature.tetrahedron.NewtonCotesOpen(0),
    quadrature.tetrahedron.NewtonCotesOpen(1),
    quadrature.tetrahedron.NewtonCotesOpen(2),
    quadrature.tetrahedron.NewtonCotesOpen(3),
    quadrature.tetrahedron.NewtonCotesOpen(4),
    quadrature.tetrahedron.NewtonCotesOpen(5),
    quadrature.tetrahedron.NewtonCotesOpen(6),
    quadrature.tetrahedron.Zienkiewicz(4),
    quadrature.tetrahedron.Zienkiewicz(5),
    quadrature.tetrahedron.ShunnHam(1),
    quadrature.tetrahedron.ShunnHam(2),
    quadrature.tetrahedron.ShunnHam(3),
    quadrature.tetrahedron.ShunnHam(4),
    quadrature.tetrahedron.ShunnHam(5),
    quadrature.tetrahedron.ShunnHam(6),
    quadrature.tetrahedron.ZhangCuiLiu(1),
    quadrature.tetrahedron.ZhangCuiLiu(2),
    quadrature.tetrahedron.Yu(1),
    quadrature.tetrahedron.Yu(2),
    quadrature.tetrahedron.Yu(3),
    quadrature.tetrahedron.Yu(4),
    quadrature.tetrahedron.Yu(5),
    quadrature.tetrahedron.HammerMarloweStroud(1),
    quadrature.tetrahedron.HammerMarloweStroud(2),
    quadrature.tetrahedron.HammerMarloweStroud(3),
    quadrature.tetrahedron.LiuVinokur(1),
    quadrature.tetrahedron.LiuVinokur(2),
    quadrature.tetrahedron.LiuVinokur(3),
    quadrature.tetrahedron.LiuVinokur(4),
    quadrature.tetrahedron.LiuVinokur(5),
    quadrature.tetrahedron.LiuVinokur(6),
    quadrature.tetrahedron.LiuVinokur(7),
    quadrature.tetrahedron.LiuVinokur(8),
    quadrature.tetrahedron.LiuVinokur(9),
    quadrature.tetrahedron.LiuVinokur(10),
    quadrature.tetrahedron.LiuVinokur(11),
    quadrature.tetrahedron.LiuVinokur(12),
    quadrature.tetrahedron.LiuVinokur(13),
    quadrature.tetrahedron.LiuVinokur(14),
    ])
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    tetrahedron = numpy.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        ])
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for k in create_monomial_exponents3(degree):
            def poly(x):
                return x[0]**k[0] * x[1]**k[1] * x[2]**k[2]
            # exact_val = _integrate_exact(poly, tetrahedron)
            exact_val = _integrate_monomial_over_standard_tet(k)
            val = quadrature.tetrahedron.integrate(
                    poly, tetrahedron, scheme
                    )
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
    tet = numpy.array([
        [numpy.cos(0.5*numpy.pi), numpy.sin(0.5*numpy.pi), -0.5],
        [numpy.cos(7.0/6.0*numpy.pi), numpy.sin(7.0/6.0*numpy.pi), -0.5],
        [numpy.cos(11.0/6.0*numpy.pi), numpy.sin(11.0/6.0*numpy.pi), -0.5],
        [0.0, 0.0, 1.0]
        ])
    quadrature.tetrahedron.show(
        tet,
        # quadrature.tetrahedron.Keast(0)
        quadrature.tetrahedron.Keast(7)
        # quadrature.tetrahedron.NewtonCotesClosed(6)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
    scheme = quadrature.tetrahedron.Keast(10)
    test_scheme(scheme)
