# -*- coding: utf-8 -*-
#
import numpy
import numpy.testing
import pytest
import quadrature
import sympy

from test_tetrahedron import _create_monomials

# import os
# import matplotlib as mpl
# if 'DISPLAY' not in os.environ:
#     # headless mode, for remote executions (and travis)
#     mpl.use('Agg')
# from matplotlib import pyplot as plt


def _integrate_exact(f, midpoint, radius):
    #
    # New coordinates: phi, theta.
    #
    phi = sympy.Symbol('phi')
    theta = sympy.Symbol('theta')
    x_xi = [
        midpoint[0] + radius*sympy.sin(phi)*sympy.cos(theta),
        midpoint[1] + radius*sympy.sin(phi)*sympy.sin(theta),
        midpoint[2] + radius*sympy.cos(phi),
        ]
    rtheta_x_rphi = sympy.sin(phi) * radius**2
    exact = sympy.integrate(
        sympy.integrate(rtheta_x_rphi * f(x_xi), (phi, 0.0, numpy.pi)),
        (theta, 0, 2*numpy.pi)
        )
    return float(exact)


@pytest.mark.parametrize('scheme', [
    quadrature.sphere.Lebedev(1),
    quadrature.sphere.Lebedev(2),
    quadrature.sphere.Lebedev(3),
    quadrature.sphere.Lebedev(4),
    quadrature.sphere.Lebedev(5),
    ])
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    midpoint = numpy.array([0.0, 0.0, 0.0])
    radius = 1.0
    success = True
    degree = 0
    max_degree = scheme.degree + 1
    while success:
        for poly in _create_monomials(degree):
            exact_val = _integrate_exact(poly, midpoint, radius)
            val = quadrature.sphere.integrate(
                    poly, midpoint, radius, scheme
                    )
            print(exact_val)
            print(val)
            if abs(exact_val - val) > 1.0e-10:
                success = False
                break
        if not success:
            break
        if degree >= max_degree:
            break
        degree += 1
    assert degree-1 >= scheme.degree
    # numpy.testing.assert_equal(degree-1, scheme.degree)
    return


# def test_show():
#     sphere = numpy.array([
#         [numpy.cos(0.5*numpy.pi), numpy.sin(0.5*numpy.pi)],
#         [numpy.cos(7.0/6.0*numpy.pi), numpy.sin(7.0/6.0*numpy.pi)],
#         [numpy.cos(11.0/6.0*numpy.pi), numpy.sin(11.0/6.0*numpy.pi)],
#         ])
#     quadrature.sphere.show(
#         sphere,
#         # quadrature.sphere.Centroid()
#         # quadrature.sphere.Vertex()
#         # quadrature.sphere.SevenPoint()
#         # quadrature.sphere.Strang(9)
#         quadrature.sphere.Dunavant(20)
#         )
#     return
#
#
if __name__ == '__main__':
    # test_show()
    # plt.show()
    scheme = quadrature.sphere.Lebedev(3)
    test_scheme(scheme)
