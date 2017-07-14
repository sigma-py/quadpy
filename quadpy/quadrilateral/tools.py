# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from .stroud import Stroud


def show(
        scheme,
        quad=numpy.array([[0, 0], [1, 0], [1, 1], [0, 1]]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given quad. The area of the disks
    around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt

    plt.plot(quad[:, 0], quad[:, 1], '-k')
    plt.plot(
        [quad[-1, 0], quad[0, 0]],
        [quad[-1, 1], quad[0, 1]],
        '-k'
        )

    plt.axis('equal')

    if not show_axes:
        plt.gca().set_axis_off()

    xi = scheme.points[:, 0]
    eta = scheme.points[:, 1]
    transformed_pts = \
        + numpy.outer(0.25 * (1.0 - xi)*(1.0 - eta), quad[0]) \
        + numpy.outer(0.25 * (1.0 + xi)*(1.0 - eta), quad[1]) \
        + numpy.outer(0.25 * (1.0 + xi)*(1.0 + eta), quad[2]) \
        + numpy.outer(0.25 * (1.0 - xi)*(1.0 + eta), quad[3])

    vol = integrate(lambda x: 1.0, quad, Stroud(1))
    helpers.plot_disks(
        plt, transformed_pts, scheme.weights, vol
        )
    plt.show()
    return


def integrate(f, quad, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points.T
    x = \
        + numpy.multiply.outer(0.25*(1.0-xi[0])*(1.0-xi[1]), quad[0]) \
        + numpy.multiply.outer(0.25*(1.0+xi[0])*(1.0-xi[1]), quad[1]) \
        + numpy.multiply.outer(0.25*(1.0+xi[0])*(1.0+xi[1]), quad[2]) \
        + numpy.multiply.outer(0.25*(1.0-xi[0])*(1.0+xi[1]), quad[3])
    x = x.T

    J0 = \
        - numpy.multiply.outer(0.25*(1-xi[1]), quad[0]) \
        + numpy.multiply.outer(0.25*(1-xi[1]), quad[1]) \
        + numpy.multiply.outer(0.25*(1+xi[1]), quad[2]) \
        - numpy.multiply.outer(0.25*(1+xi[1]), quad[3])
    J0 = J0.T
    J1 = \
        - numpy.multiply.outer(0.25*(1-xi[0]), quad[0]) \
        - numpy.multiply.outer(0.25*(1+xi[0]), quad[1]) \
        + numpy.multiply.outer(0.25*(1+xi[0]), quad[2]) \
        + numpy.multiply.outer(0.25*(1-xi[0]), quad[3])
    J1 = J1.T
    det = (J0[0]*J1[1] - J1[0]*J0[1])

    return sumfun(scheme.weights * f(x) * abs(det), axis=-1)
