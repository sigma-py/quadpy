# -*- coding: utf-8 -*-
#
import numpy
import numpy.testing
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


def _create_monomials(degree):
    '''Returns a list of all monomials of degree :degree:.
    '''
    return [
        lambda x: x[0]**(degree-i-j) * x[1]**i * x[2]**j
        for i in range(degree+1)
        for j in range(degree-i+1)
        ]


def test_generator():
    tetrahedron = numpy.array([
        [-1.0, -2.0, -0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 2.0],
        ])
    schemes = [
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
        quadrature.tetrahedron.NewtonCotesClosed(0),
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
        ]
    for scheme in schemes:
        yield check_scheme, scheme, tetrahedron


def check_scheme(scheme, tetrahedron):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    success = True
    degree = 0
    max_degree = 100
    while success:
        for poly in _create_monomials(degree):
            exact_val = _integrate_exact(poly, tetrahedron)
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
        # quadrature.tetrahedron.Keast(7)
        quadrature.tetrahedron.NewtonCotesClosed(6)
        )
    return


if __name__ == '__main__':
    test_show()
    plt.show()
