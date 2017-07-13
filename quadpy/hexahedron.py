# -*- coding: utf-8 -*-
#
import numpy

from . import helpers


def show(
        scheme,
        hexa=numpy.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
            ]),
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
        [hexa[0], hexa[1]],
        [hexa[1], hexa[2]],
        [hexa[2], hexa[3]],
        [hexa[3], hexa[0]],
        #
        [hexa[4], hexa[5]],
        [hexa[5], hexa[6]],
        [hexa[6], hexa[7]],
        [hexa[7], hexa[4]],
        #
        [hexa[0], hexa[4]],
        [hexa[1], hexa[5]],
        [hexa[2], hexa[6]],
        [hexa[3], hexa[7]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    xi = scheme.points[:, 0]
    eta = scheme.points[:, 1]
    zeta = scheme.points[:, 2]
    transformed_pts = \
        + numpy.outer(0.125 * (1.0 - xi)*(1.0 - eta)*(1.0 - zeta), hexa[0]) \
        + numpy.outer(0.125 * (1.0 + xi)*(1.0 - eta)*(1.0 - zeta), hexa[1]) \
        + numpy.outer(0.125 * (1.0 + xi)*(1.0 + eta)*(1.0 - zeta), hexa[2]) \
        + numpy.outer(0.125 * (1.0 - xi)*(1.0 + eta)*(1.0 - zeta), hexa[3]) \
        + numpy.outer(0.125 * (1.0 - xi)*(1.0 - eta)*(1.0 + zeta), hexa[4]) \
        + numpy.outer(0.125 * (1.0 + xi)*(1.0 - eta)*(1.0 + zeta), hexa[5]) \
        + numpy.outer(0.125 * (1.0 + xi)*(1.0 + eta)*(1.0 + zeta), hexa[6]) \
        + numpy.outer(0.125 * (1.0 - xi)*(1.0 + eta)*(1.0 + zeta), hexa[7])

    vol = integrate(lambda x: 1.0, hexa, scheme)
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol
        )

    plt.show()
    return


def integrate(f, hexa, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points.T
    mo = numpy.multiply.outer
    x = \
        + mo(0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2]), hexa[0]) \
        + mo(0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2]), hexa[1]) \
        + mo(0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2]), hexa[2]) \
        + mo(0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2]), hexa[3]) \
        + mo(0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2]), hexa[4]) \
        + mo(0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2]), hexa[5]) \
        + mo(0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2]), hexa[6]) \
        + mo(0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2]), hexa[7])
    x = x.T

    J0 = \
        - numpy.multiply.outer(0.125*(1 - xi[1])*(1.0 - xi[2]), hexa[0]) \
        + numpy.multiply.outer(0.125*(1 - xi[1])*(1.0 - xi[2]), hexa[1]) \
        + numpy.multiply.outer(0.125*(1 + xi[1])*(1.0 - xi[2]), hexa[2]) \
        - numpy.multiply.outer(0.125*(1 + xi[1])*(1.0 - xi[2]), hexa[3]) \
        - numpy.multiply.outer(0.125*(1 - xi[1])*(1.0 + xi[2]), hexa[4]) \
        + numpy.multiply.outer(0.125*(1 - xi[1])*(1.0 + xi[2]), hexa[5]) \
        + numpy.multiply.outer(0.125*(1 + xi[1])*(1.0 + xi[2]), hexa[6]) \
        - numpy.multiply.outer(0.125*(1 + xi[1])*(1.0 + xi[2]), hexa[7])
    J0 = J0.T
    J1 = \
        - numpy.multiply.outer(0.125*(1 - xi[0])*(1.0 - xi[2]), hexa[0]) \
        - numpy.multiply.outer(0.125*(1 + xi[0])*(1.0 - xi[2]), hexa[1]) \
        + numpy.multiply.outer(0.125*(1 + xi[0])*(1.0 - xi[2]), hexa[2]) \
        + numpy.multiply.outer(0.125*(1 - xi[0])*(1.0 - xi[2]), hexa[3]) \
        - numpy.multiply.outer(0.125*(1 - xi[0])*(1.0 + xi[2]), hexa[4]) \
        - numpy.multiply.outer(0.125*(1 + xi[0])*(1.0 + xi[2]), hexa[5]) \
        + numpy.multiply.outer(0.125*(1 + xi[0])*(1.0 + xi[2]), hexa[6]) \
        + numpy.multiply.outer(0.125*(1 - xi[0])*(1.0 + xi[2]), hexa[7])
    J1 = J1.T
    J2 = \
        - numpy.multiply.outer(0.125*(1.0 - xi[0])*(1.0 - xi[1]), hexa[0]) \
        - numpy.multiply.outer(0.125*(1.0 + xi[0])*(1.0 - xi[1]), hexa[1]) \
        - numpy.multiply.outer(0.125*(1.0 + xi[0])*(1.0 + xi[1]), hexa[2]) \
        - numpy.multiply.outer(0.125*(1.0 - xi[0])*(1.0 + xi[1]), hexa[3]) \
        + numpy.multiply.outer(0.125*(1.0 - xi[0])*(1.0 - xi[1]), hexa[4]) \
        + numpy.multiply.outer(0.125*(1.0 + xi[0])*(1.0 - xi[1]), hexa[5]) \
        + numpy.multiply.outer(0.125*(1.0 + xi[0])*(1.0 + xi[1]), hexa[6]) \
        + numpy.multiply.outer(0.125*(1.0 - xi[0])*(1.0 + xi[1]), hexa[7])
    J2 = J2.T
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]

    return sumfun(scheme.weights * f(x) * abs(det), axis=-1)


class Product(object):
    def __init__(self, scheme1d):
        self.schemes = \
            scheme1d if isinstance(scheme1d, list) \
            else 3 * [scheme1d]

        wy, wz, wx = numpy.meshgrid(
            self.schemes[0].weights,
            self.schemes[1].weights,
            self.schemes[2].weights,
            )
        weights = numpy.vstack([
            wx.flatten(),
            wy.flatten(),
            wz.flatten()
            ]).T
        self.weights = numpy.prod(weights, axis=1)
        # the order, yeah...
        y, z, x = numpy.meshgrid(
            self.schemes[0].points,
            self.schemes[1].points,
            self.schemes[2].points,
            )
        self.points = numpy.vstack([
            x.flatten(),
            y.flatten(),
            z.flatten(),
            ]).T

        self.degree = min([s.degree for s in self.schemes])
        return
