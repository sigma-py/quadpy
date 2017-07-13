# -*- coding: utf-8 -*-
#
import numpy

from .felippa import Felippa
from .. import helpers


def show(
        scheme,
        pyra=numpy.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0.5, 0.5, 1.0],
            ])
        ):
    '''Shows the quadrature points on a given pyramid. The size of the
    balls around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    edges = numpy.array([
        [pyra[0], pyra[1]],
        [pyra[1], pyra[2]],
        [pyra[2], pyra[3]],
        [pyra[3], pyra[0]],
        #
        [pyra[0], pyra[4]],
        [pyra[1], pyra[4]],
        [pyra[2], pyra[4]],
        [pyra[3], pyra[4]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    xi = scheme.points[:, 0]
    eta = scheme.points[:, 1]
    zeta = scheme.points[:, 2]
    transformed_pts = \
        + numpy.outer(pyra[0], 0.125*(1.0-xi)*(1.0-eta)*(1-zeta)) \
        + numpy.outer(pyra[1], 0.125*(1.0+xi)*(1.0-eta)*(1-zeta)) \
        + numpy.outer(pyra[2], 0.125*(1.0+xi)*(1.0+eta)*(1-zeta)) \
        + numpy.outer(pyra[3], 0.125*(1.0-xi)*(1.0+eta)*(1-zeta)) \
        + numpy.outer(pyra[4], 0.500*(1.0+zeta))
    transformed_pts = transformed_pts.T

    vol = integrate(lambda x: 1.0, pyra, Felippa(1))
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol
        )
    plt.show()
    return


def _get_det_J(pyra, xi):
    J0 = \
        - numpy.multiply.outer(0.125*(1.0-xi[1])*(1-xi[2]), pyra[0]) \
        + numpy.multiply.outer(0.125*(1.0-xi[1])*(1-xi[2]), pyra[1]) \
        + numpy.multiply.outer(0.125*(1.0+xi[1])*(1-xi[2]), pyra[2]) \
        - numpy.multiply.outer(0.125*(1.0+xi[1])*(1-xi[2]), pyra[3])
    J0 = J0.T
    J1 = \
        - numpy.multiply.outer(0.125*(1.0-xi[0])*(1-xi[2]), pyra[0]) \
        - numpy.multiply.outer(0.125*(1.0+xi[0])*(1-xi[2]), pyra[1]) \
        + numpy.multiply.outer(0.125*(1.0+xi[0])*(1-xi[2]), pyra[2]) \
        + numpy.multiply.outer(0.125*(1.0-xi[0])*(1-xi[2]), pyra[3])
    J1 = J1.T
    J2 = \
        - numpy.multiply.outer(0.125*(1.0-xi[0])*(1.0-xi[1]), pyra[0]) \
        - numpy.multiply.outer(0.125*(1.0+xi[0])*(1.0-xi[1]), pyra[1]) \
        - numpy.multiply.outer(0.125*(1.0+xi[0])*(1.0+xi[1]), pyra[2]) \
        - numpy.multiply.outer(0.125*(1.0-xi[0])*(1.0+xi[1]), pyra[3]) \
        + numpy.multiply.outer(0.500*numpy.ones(1), pyra[4])
    J2 = J2.T
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]
    return det.T


def integrate(f, pyra, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points.T
    mo = numpy.multiply.outer
    x = \
        + mo(0.125*(1.0-xi[0])*(1.0-xi[1])*(1-xi[2]), pyra[0]) \
        + mo(0.125*(1.0+xi[0])*(1.0-xi[1])*(1-xi[2]), pyra[1]) \
        + mo(0.125*(1.0+xi[0])*(1.0+xi[1])*(1-xi[2]), pyra[2]) \
        + mo(0.125*(1.0-xi[0])*(1.0+xi[1])*(1-xi[2]), pyra[3]) \
        + mo(0.500*(1.0+xi[2]), pyra[4])
    x = x.T
    det = _get_det_J(pyra, xi)
    return sumfun(scheme.weights * f(x) * abs(det.T), axis=-1)
