# -*- coding: utf-8 -*-
#
import math
import numpy
from . import tetrahedron


def volume(hexa):
    # Split the hexahedron into five tetrahedra, cf.
    # <http://www.cvel.clemson.edu/modeling/EMAG/EMAP/emap4/Meshing.html>,
    # and get the volume of those. The numbering is assumed to be
    #
    #     7___6
    #    /    /|
    #  4/___5/ |
    #   |    | /2
    #   |____|/
    #  0      1
    #
    # (3 is on the back side below 7).
    tets = numpy.array([
        [hexa[0], hexa[1], hexa[2], hexa[5]],
        [hexa[1], hexa[2], hexa[3], hexa[7]],
        [hexa[0], hexa[4], hexa[5], hexa[7]],
        [hexa[2], hexa[5], hexa[6], hexa[7]],
        [hexa[0], hexa[2], hexa[5], hexa[7]],
        ])
    vol = numpy.sum([tetrahedron.volume(tet) for tet in tets])
    return vol


def show(hexa, scheme, ball_scale=1.0, alpha=0.3):
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

    hexa_vol = volume(hexa)
    phi, theta = numpy.mgrid[0:numpy.pi:101j, 0:2*numpy.pi:101j]
    x = numpy.sin(phi)*numpy.cos(theta)
    y = numpy.sin(phi)*numpy.sin(theta)
    z = numpy.cos(phi)
    for tp, weight in zip(transformed_pts, scheme.weights):
        color = 'b' if weight >= 0 else 'r'
        # highlight ball center
        plt.plot([tp[0]], [tp[1]], [tp[2]], '.' + color)

        # plot ball
        # scale the circle volume according to the weight
        ref_vol = 8.0
        r = ball_scale * (
            hexa_vol * abs(weight) / ref_vol / (4.0/3.0 * numpy.pi)
            )**(1.0/3.0)

        ax.plot_surface(
            r*x + tp[0], r*y + tp[1], r*z + tp[2],
            color=color,
            alpha=alpha,
            linewidth=0
            )

    # http://stackoverflow.com/a/21765085/353337
    alpha = 1.3
    max_range = alpha * 0.5 * numpy.array([
        hexa[:, 0].max() - hexa[:, 0].min(),
        hexa[:, 1].max() - hexa[:, 1].min(),
        hexa[:, 2].max() - hexa[:, 2].min(),
        ]).max()
    mid_x = 0.5 * (hexa[:, 0].max() + hexa[:, 0].min())
    mid_y = 0.5 * (hexa[:, 1].max() + hexa[:, 1].min())
    mid_z = 0.5 * (hexa[:, 2].max() + hexa[:, 2].min())
    #
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

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
