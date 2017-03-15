# -*- coding: utf-8 -*-
#
import math
import numpy

from . import helpers


def show(hexa, scheme):
    '''Shows the quadrature points on a given hexahedron. The size of the
    balls around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

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

    vol = integrate(lambda x: numpy.ones(1), hexa, scheme)
    helpers.plot_balls(
        plt, ax, transformed_pts, scheme.weights, vol,
        hexa[:, 0].min(), hexa[:, 0].max(),
        hexa[:, 1].min(), hexa[:, 1].max(),
        hexa[:, 2].min(), hexa[:, 2].max(),
        )

    plt.show()
    return


def integrate(f, hexa, scheme):
    xi = scheme.points.T
    x = \
        + numpy.outer(hexa[0], 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0-xi[2])) \
        + numpy.outer(hexa[1], 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0-xi[2])) \
        + numpy.outer(hexa[2], 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0-xi[2])) \
        + numpy.outer(hexa[3], 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0-xi[2])) \
        + numpy.outer(hexa[4], 0.125*(1.0-xi[0])*(1.0-xi[1])*(1.0+xi[2])) \
        + numpy.outer(hexa[5], 0.125*(1.0+xi[0])*(1.0-xi[1])*(1.0+xi[2])) \
        + numpy.outer(hexa[6], 0.125*(1.0+xi[0])*(1.0+xi[1])*(1.0+xi[2])) \
        + numpy.outer(hexa[7], 0.125*(1.0-xi[0])*(1.0+xi[1])*(1.0+xi[2]))
    J0 = \
        - numpy.outer(hexa[0], 0.125*(1 - xi[1])*(1.0 - xi[2])) \
        + numpy.outer(hexa[1], 0.125*(1 - xi[1])*(1.0 - xi[2])) \
        + numpy.outer(hexa[2], 0.125*(1 + xi[1])*(1.0 - xi[2])) \
        - numpy.outer(hexa[3], 0.125*(1 + xi[1])*(1.0 - xi[2])) \
        - numpy.outer(hexa[4], 0.125*(1 - xi[1])*(1.0 + xi[2])) \
        + numpy.outer(hexa[5], 0.125*(1 - xi[1])*(1.0 + xi[2])) \
        + numpy.outer(hexa[6], 0.125*(1 + xi[1])*(1.0 + xi[2])) \
        - numpy.outer(hexa[7], 0.125*(1 + xi[1])*(1.0 + xi[2]))
    J1 = \
        - numpy.outer(hexa[0], 0.125*(1 - xi[0])*(1.0 - xi[2])) \
        - numpy.outer(hexa[1], 0.125*(1 + xi[0])*(1.0 - xi[2])) \
        + numpy.outer(hexa[2], 0.125*(1 + xi[0])*(1.0 - xi[2])) \
        + numpy.outer(hexa[3], 0.125*(1 - xi[0])*(1.0 - xi[2])) \
        - numpy.outer(hexa[4], 0.125*(1 - xi[0])*(1.0 + xi[2])) \
        - numpy.outer(hexa[5], 0.125*(1 + xi[0])*(1.0 + xi[2])) \
        + numpy.outer(hexa[6], 0.125*(1 + xi[0])*(1.0 + xi[2])) \
        + numpy.outer(hexa[7], 0.125*(1 - xi[0])*(1.0 + xi[2]))
    J2 = \
        - numpy.outer(hexa[0], 0.125*(1.0 - xi[0])*(1.0 - xi[1])) \
        - numpy.outer(hexa[1], 0.125*(1.0 + xi[0])*(1.0 - xi[1])) \
        - numpy.outer(hexa[2], 0.125*(1.0 + xi[0])*(1.0 + xi[1])) \
        - numpy.outer(hexa[3], 0.125*(1.0 - xi[0])*(1.0 + xi[1])) \
        + numpy.outer(hexa[4], 0.125*(1.0 - xi[0])*(1.0 - xi[1])) \
        + numpy.outer(hexa[5], 0.125*(1.0 + xi[0])*(1.0 - xi[1])) \
        + numpy.outer(hexa[6], 0.125*(1.0 + xi[0])*(1.0 + xi[1])) \
        + numpy.outer(hexa[7], 0.125*(1.0 - xi[0])*(1.0 + xi[1]))
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]

    return math.fsum(scheme.weights * f(x).T * abs(det))


class From1d(object):
    def __init__(self, scheme1d):
        wy, wz, wx = numpy.meshgrid(
            scheme1d.weights, scheme1d.weights, scheme1d.weights
            )
        weights = numpy.vstack([
            wx.flatten(),
            wy.flatten(),
            wz.flatten()
            ]).T
        self.weights = numpy.prod(weights, axis=1)
        # the order, yeah...
        y, z, x = numpy.meshgrid(
            scheme1d.points, scheme1d.points, scheme1d.points
            )
        self.points = numpy.vstack([
            x.flatten(),
            y.flatten(),
            z.flatten(),
            ]).T

        self.degree = scheme1d.degree
        return
