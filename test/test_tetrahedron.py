# -*- coding: utf-8 -*-
#
from helpers import create_monomial_exponents3, check_degree
import math
from matplotlib import pyplot as plt
import numpy
import pytest
import quadpy
import sympy


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
    abs_det_J = 6 * quadpy.tetrahedron.volume(tetrahedron)
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

    \int_T x_0^k0 * x1^k1 * x2^k2 = (k0!*k1!*k2!) / (4+k0+k1+k2)!,

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


@pytest.mark.parametrize(
    'scheme',
    [quadpy.tetrahedron.HammerMarloweStroud(k) for k in [1, 2, 3]]
    + [quadpy.tetrahedron.NewtonCotesClosed(k) for k in range(1, 7)]
    + [quadpy.tetrahedron.NewtonCotesOpen(k) for k in range(7)]
    + [quadpy.tetrahedron.Yu(k) for k in range(1, 6)]
    + [quadpy.tetrahedron.Keast(k) for k in range(11)]
    + [quadpy.tetrahedron.LiuVinokur(k) for k in range(1, 15)]
    + [quadpy.tetrahedron.Zienkiewicz(k) for k in [4, 5]]
    + [quadpy.tetrahedron.ZhangCuiLiu(k) for k in [1, 2]]
    + [quadpy.tetrahedron.ShunnHam(k) for k in range(1, 7)]
    + [quadpy.tetrahedron.WilliamsShunnJameson()]
    + [quadpy.tetrahedron.XiaoGimbutas(k) for k in range(1, 16)]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    tetrahedron = numpy.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        ])
    degree = check_degree(
            lambda poly: quadpy.tetrahedron.integrate(
                poly, tetrahedron, scheme
                ),
            _integrate_monomial_over_standard_tet,
            create_monomial_exponents3,
            scheme.degree + 1
            )
    assert degree == scheme.degree
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.tetrahedron.HammerMarloweStroud(3)]
    )
def test_show(scheme):
    tet = numpy.array([
        [numpy.cos(0.5*numpy.pi), numpy.sin(0.5*numpy.pi), -0.5],
        [numpy.cos(7.0/6.0*numpy.pi), numpy.sin(7.0/6.0*numpy.pi), -0.5],
        [numpy.cos(11.0/6.0*numpy.pi), numpy.sin(11.0/6.0*numpy.pi), -0.5],
        [0.0, 0.0, 1.0]
        ])
    quadpy.tetrahedron.show(scheme, tet)
    return


if __name__ == '__main__':
    scheme = quadpy.tetrahedron.ShunnHam(2)
    test_scheme(scheme)
    test_show(scheme)
    plt.show()
