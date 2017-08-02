# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from . import felippa
from .. import helpers


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(
        scheme,
        wedge=numpy.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [0, 1, 1],
            ]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given wedge. The size of the
    balls around the points coincides with their weights.
    '''
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    if not show_axes:
        ax.set_axis_off()

    edges = numpy.array([
        [wedge[0], wedge[1]],
        [wedge[1], wedge[2]],
        [wedge[0], wedge[2]],
        #
        [wedge[3], wedge[4]],
        [wedge[4], wedge[5]],
        [wedge[5], wedge[3]],
        #
        [wedge[0], wedge[3]],
        [wedge[1], wedge[4]],
        [wedge[2], wedge[5]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    xi = scheme.points[:, 0]
    eta = scheme.points[:, 1]
    zeta = scheme.points[:, 2]
    transformed_pts = \
        + numpy.outer(0.5 * (1.0 - xi - eta)*(1.0 - zeta), wedge[0]) \
        + numpy.outer(0.5 * xi * (1.0 - zeta), wedge[1]) \
        + numpy.outer(0.5 * eta * (1.0 - zeta), wedge[2]) \
        + numpy.outer(0.5 * (1.0 - xi - eta)*(1.0 - zeta), wedge[3]) \
        + numpy.outer(0.5 * xi * (1.0 + zeta), wedge[4]) \
        + numpy.outer(0.5 * eta * (1.0 + zeta), wedge[5])

    vol = integrate(lambda x: numpy.ones(1), wedge, felippa.Felippa(1))
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol
        )
    plt.show()
    return


def integrate(f, wedge, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points.T
    mo = numpy.multiply.outer
    x = \
        + mo(0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2]), wedge[0]) \
        + mo(0.5 * xi[0] * (1.0-xi[2]), wedge[1]) \
        + mo(0.5 * xi[1] * (1.0-xi[2]), wedge[2]) \
        + mo(0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2]), wedge[3]) \
        + mo(0.5 * xi[0] * (1.0+xi[2]), wedge[4]) \
        + mo(0.5 * xi[1] * (1.0+xi[2]), wedge[5])
    x = x.T
    J0 = \
        - mo(0.5*(1.0 - xi[2]), wedge[0]) \
        + mo(0.5*(1.0 - xi[2]), wedge[2]) \
        - mo(0.5*(1.0 + xi[2]), wedge[3]) \
        + mo(0.5*(1.0 + xi[2]), wedge[5])
    J0 = J0.T
    J1 = \
        - mo(0.5*(1.0 - xi[2]), wedge[0]) \
        + mo(0.5*(1.0 - xi[2]), wedge[1]) \
        - mo(0.5*(1.0 + xi[2]), wedge[3]) \
        + mo(0.5*(1.0 + xi[2]), wedge[4])
    J1 = J1.T
    J2 = \
        - mo(0.5 * (1.0-xi[0]-xi[1]), wedge[0]) \
        - mo(0.5 * xi[0], wedge[1]) \
        - mo(0.5 * xi[1], wedge[2]) \
        + mo(0.5 * (1.0-xi[0]-xi[1]), wedge[3]) \
        + mo(0.5 * xi[0], wedge[4]) \
        + mo(0.5 * xi[1], wedge[5])
    J2 = J2.T
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]

    return sumfun(scheme.weights * f(x) * abs(det), axis=-1)
