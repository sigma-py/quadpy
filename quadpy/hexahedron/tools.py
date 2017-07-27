# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def cube_points(x, y, z):
    '''Given the end points of a cube aligned with the coordinate axes, this
    returns the corner points of the cube in the correct data structure.
    '''
    return numpy.moveaxis(
            numpy.array(numpy.meshgrid(x, y, z, indexing='ij')),
            0, -1
            )


def show(
        scheme,
        hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given hexahedron. The size of the
    balls around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    if not show_axes:
        ax.set_axis_off()

    edges = numpy.array([
        [hexa[0][0][0], hexa[1, 0, 0]],
        [hexa[1][0][0], hexa[1, 1, 0]],
        [hexa[1][1][0], hexa[0, 1, 0]],
        [hexa[0][1][0], hexa[0, 0, 0]],
        #
        [hexa[0][0][1], hexa[1, 0, 1]],
        [hexa[1][0][1], hexa[1, 1, 1]],
        [hexa[1][1][1], hexa[0, 1, 1]],
        [hexa[0][1][1], hexa[0, 0, 1]],
        #
        [hexa[0][0][0], hexa[0, 0, 1]],
        [hexa[1][0][0], hexa[1, 0, 1]],
        [hexa[1][1][0], hexa[1, 1, 1]],
        [hexa[0][1][0], hexa[0, 1, 1]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    transformed_pts = _transform_to_hexa(scheme.points.T, hexa)

    vol = integrate(lambda x: 1.0, hexa, scheme)
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol
        )

    plt.show()
    return


def _transform_to_hexa(xi, hexa):
    '''Transform the points xi from the reference hexahedron to the hexahedron
    `hexa`.
    '''
    mo = numpy.multiply.outer
    return (
        + mo(0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2]), hexa[0, 0, 0])
        + mo(0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2]), hexa[1, 0, 0])
        + mo(0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2]), hexa[1, 1, 0])
        + mo(0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2]), hexa[0, 1, 0])
        + mo(0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2]), hexa[0, 0, 1])
        + mo(0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2]), hexa[1, 0, 1])
        + mo(0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2]), hexa[1, 1, 1])
        + mo(0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2]), hexa[0, 1, 1])
        )
    return


def _get_detJ(xi, hexa):
    '''Get the determinant of the transformation matrix.
    '''
    J0 = \
        - numpy.multiply.outer(0.125*(1 - xi[1])*(1 - xi[2]), hexa[0, 0, 0]) \
        + numpy.multiply.outer(0.125*(1 - xi[1])*(1 - xi[2]), hexa[1, 0, 0]) \
        + numpy.multiply.outer(0.125*(1 + xi[1])*(1 - xi[2]), hexa[1, 1, 0]) \
        - numpy.multiply.outer(0.125*(1 + xi[1])*(1 - xi[2]), hexa[0, 1, 0]) \
        - numpy.multiply.outer(0.125*(1 - xi[1])*(1 + xi[2]), hexa[0, 0, 1]) \
        + numpy.multiply.outer(0.125*(1 - xi[1])*(1 + xi[2]), hexa[1, 0, 1]) \
        + numpy.multiply.outer(0.125*(1 + xi[1])*(1 + xi[2]), hexa[1, 1, 1]) \
        - numpy.multiply.outer(0.125*(1 + xi[1])*(1 + xi[2]), hexa[0, 1, 1])
    J0 = J0.T
    J1 = \
        - numpy.multiply.outer(0.125*(1 - xi[0])*(1 - xi[2]), hexa[0, 0, 0]) \
        - numpy.multiply.outer(0.125*(1 + xi[0])*(1 - xi[2]), hexa[1, 0, 0]) \
        + numpy.multiply.outer(0.125*(1 + xi[0])*(1 - xi[2]), hexa[1, 1, 0]) \
        + numpy.multiply.outer(0.125*(1 - xi[0])*(1 - xi[2]), hexa[0, 1, 0]) \
        - numpy.multiply.outer(0.125*(1 - xi[0])*(1 + xi[2]), hexa[0, 0, 1]) \
        - numpy.multiply.outer(0.125*(1 + xi[0])*(1 + xi[2]), hexa[1, 0, 1]) \
        + numpy.multiply.outer(0.125*(1 + xi[0])*(1 + xi[2]), hexa[1, 1, 1]) \
        + numpy.multiply.outer(0.125*(1 - xi[0])*(1 + xi[2]), hexa[0, 1, 1])
    J1 = J1.T
    J2 = \
        - numpy.multiply.outer(0.125*(1 - xi[0])*(1 - xi[1]), hexa[0, 0, 0]) \
        - numpy.multiply.outer(0.125*(1 + xi[0])*(1 - xi[1]), hexa[1, 0, 0]) \
        - numpy.multiply.outer(0.125*(1 + xi[0])*(1 + xi[1]), hexa[1, 1, 0]) \
        - numpy.multiply.outer(0.125*(1 - xi[0])*(1 + xi[1]), hexa[0, 1, 0]) \
        + numpy.multiply.outer(0.125*(1 - xi[0])*(1 - xi[1]), hexa[0, 0, 1]) \
        + numpy.multiply.outer(0.125*(1 + xi[0])*(1 - xi[1]), hexa[1, 0, 1]) \
        + numpy.multiply.outer(0.125*(1 + xi[0])*(1 + xi[1]), hexa[1, 1, 1]) \
        + numpy.multiply.outer(0.125*(1 - xi[0])*(1 + xi[1]), hexa[0, 1, 1])
    J2 = J2.T
    return (
        + J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2]
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]
        )


def integrate(f, hexa, scheme, sumfun=helpers.kahan_sum):
    x = _transform_to_hexa(scheme.points.T, hexa).T
    detJ = _get_detJ(scheme.points.T, hexa)
    return sumfun(scheme.weights * f(x) * abs(detJ), axis=-1)
