# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from . import felippa
from .. import helpers


def integrate(f, wedge, scheme, sumfun=helpers.kahan_sum):
    x = _transform(scheme.points.T, wedge)
    det = _get_detJ(scheme.points.T, wedge)
    return sumfun(scheme.weights * f(x) * abs(det), axis=-1)


def _transform(xi, wedge):
    mo = numpy.multiply.outer
    return (
        + mo(0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2]), wedge[0])
        + mo(0.5 * xi[0] * (1.0-xi[2]), wedge[1])
        + mo(0.5 * xi[1] * (1.0-xi[2]), wedge[2])
        + mo(0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2]), wedge[3])
        + mo(0.5 * xi[0] * (1.0+xi[2]), wedge[4])
        + mo(0.5 * xi[1] * (1.0+xi[2]), wedge[5])
        ).T


def _get_detJ(xi, wedge):
    mo = numpy.multiply.outer
    J0 = (
        - mo(0.5*(1.0 - xi[2]), wedge[0])
        + mo(0.5*(1.0 - xi[2]), wedge[2])
        - mo(0.5*(1.0 + xi[2]), wedge[3])
        + mo(0.5*(1.0 + xi[2]), wedge[5])
        ).T
    J1 = (
        - mo(0.5*(1.0 - xi[2]), wedge[0])
        + mo(0.5*(1.0 - xi[2]), wedge[1])
        - mo(0.5*(1.0 + xi[2]), wedge[3])
        + mo(0.5*(1.0 + xi[2]), wedge[4])
        ).T
    J2 = (
        - mo(0.5 * (1.0-xi[0]-xi[1]), wedge[0])
        - mo(0.5 * xi[0], wedge[1])
        - mo(0.5 * xi[1], wedge[2])
        + mo(0.5 * (1.0-xi[0]-xi[1]), wedge[3])
        + mo(0.5 * xi[0], wedge[4])
        + mo(0.5 * xi[1], wedge[5])
        ).T
    det = (
        + J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2]
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]
        )
    return det


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(
        scheme,
        wedge=numpy.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [0, 1, 1],
            ], dtype=float),
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
    edges = numpy.moveaxis(edges, 1, 2)
    for edge in edges:
        plt.plot(*edge, color='k')

    transformed_pts = _transform(scheme.points.T, wedge).T

    vol = integrate(lambda x: numpy.ones(1), wedge, felippa.Felippa(1))
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol
        )
    return


def show_mayavi(
        scheme,
        wedge=numpy.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [0, 1, 1],
            ], dtype=float),
        ):
    '''Shows the quadrature points on a given tetrahedron. The size of the
    balls around the points coincides with their weights.
    '''
    import mayavi.mlab as mlab

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
        mlab.plot3d(edge[:, 0], edge[:, 1], edge[:, 2], tube_radius=1.0e-2)

    helpers.plot_spheres_mayavi(
            _transform(scheme.points.T, wedge).T,
            scheme.weights,
            integrate(lambda x: numpy.ones(1), wedge, felippa.Felippa(1))
            )
    mlab.show()
    return
