# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from ..simplex import transform, get_vol


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
    # pylint: disable=relative-import, unused-variable
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

    transformed_pts = transform(scheme.points.T, tet.T)

    vol = get_vol(tet)
    helpers.plot_spheres(plt, ax, transformed_pts, scheme.weights, vol)
    plt.show()
    return
