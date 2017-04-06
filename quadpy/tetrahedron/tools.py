# -*- coding: utf-8 -*-
#
from .keast import Keast

from .. import helpers

import numpy


def show(
        scheme,
        tet=numpy.array([
            [+1, 0, -1.0/numpy.sqrt(2.0)],
            [-1, 0, -1.0/numpy.sqrt(2.0)],
            [0, +1, +1.0/numpy.sqrt(2.0)],
            [0, -1, +1.0/numpy.sqrt(2.0)],
            ]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given tetrahedron. The size of the
    balls around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    if not show_axes:
        plt.gca().set_axis_off()

    edges = numpy.array([
        [tet[0], tet[1]],
        [tet[0], tet[2]],
        [tet[0], tet[3]],
        [tet[1], tet[2]],
        [tet[1], tet[3]],
        [tet[2], tet[3]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    transformed_pts = \
        + numpy.outer(
            (1.0 - scheme.points[:, 0]
                 - scheme.points[:, 1]
                 - scheme.points[:, 2]),
            tet[0]
            ) \
        + numpy.outer(scheme.points[:, 0], tet[1]) \
        + numpy.outer(scheme.points[:, 1], tet[2]) \
        + numpy.outer(scheme.points[:, 2], tet[3])

    vol = integrate(lambda x: numpy.ones(1), tet, Keast(0))
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol,
        tet[:, 0].min(), tet[:, 0].max(),
        tet[:, 1].min(), tet[:, 1].max(),
        tet[:, 2].min(), tet[:, 2].max(),
        )
    plt.show()
    return


def integrate(f, tetrahedron, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points.T
    x = \
        + numpy.multiply.outer(1.0 - xi[0] - xi[1] - xi[2], tetrahedron[0]) \
        + numpy.multiply.outer(xi[0], tetrahedron[1]) \
        + numpy.multiply.outer(xi[1], tetrahedron[2]) \
        + numpy.multiply.outer(xi[2], tetrahedron[3])
    x = x.T

    # det is the signed volume of the tetrahedron
    J0 = (tetrahedron[1] - tetrahedron[0]).T
    J1 = (tetrahedron[2] - tetrahedron[0]).T
    J2 = (tetrahedron[3] - tetrahedron[0]).T
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]
    # reference volume
    det *= 1.0/6.0

    return sumfun(numpy.rollaxis(scheme.weights * f(x), -1) * abs(det))
