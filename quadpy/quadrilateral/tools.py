# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from .stroud import Stroud
# from ..ncube import transform


def rectangle_points(x, y):
    '''Given the end points of a rectangle aligned with the coordinate axes,
    this returns the corner points of the cube in the correct data structure.
    '''
    return numpy.moveaxis(
            numpy.array(numpy.meshgrid(x, y, indexing='ij')),
            0, -1
            )


def show(
        scheme,
        quad=rectangle_points([0.0, 1.0], [0.0, 1.0]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given quad. The area of the disks
    around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt

    plt.plot(quad[0][0], quad[1][0], '-k')
    plt.plot(quad[1][0], quad[1][1], '-k')
    plt.plot(quad[1][1], quad[0][1], '-k')
    plt.plot(quad[0][1], quad[0][0], '-k')

    plt.axis('equal')

    if not show_axes:
        plt.gca().set_axis_off()

    transformed_pts = transform(scheme.points.T, quad)

    vol = integrate(lambda x: 1.0, quad, Stroud(1))
    helpers.plot_disks(
        plt, transformed_pts, scheme.weights, vol
        )
    plt.show()
    return


def transform(xi, quad):
    '''Transform the points xi from the reference quad to the quad `quad`.
    '''
    # For d==2, the result used to be computed with
    #
    # out = (
    #     + outer(0.25*(1.0-xi[0])*(1.0-xi[1]), quad[0, 0])
    #     + outer(0.25*(1.0+xi[0])*(1.0-xi[1]), quad[1, 0])
    #     + outer(0.25*(1.0-xi[0])*(1.0+xi[1]), quad[0, 1])
    #     + outer(0.25*(1.0+xi[0])*(1.0+xi[1]), quad[1, 1])
    #     )
    #
    # This array of multiplications and additions is reminiscent of dot(), and
    # indeed tensordot() can handle the situation. We just need to compute the
    # `1+-xi` products and align them with `quad`.
    one_mp_xi = numpy.stack([
        0.5 * (1.0 - xi),
        0.5 * (1.0 + xi),
        ], axis=1)
    a = helpers.n_outer(one_mp_xi)

    # TODO kahan tensordot
    # <https://stackoverflow.com/q/45372098/353337>
    d = 2
    return numpy.tensordot(a, quad, axes=(range(d), range(d)))


def _get_detJ(xi, quad):
    '''Get the determinant of the transformation matrix.
    '''
    J0 = (
        - numpy.multiply.outer(0.25*(1-xi[1]), quad[0, 0])
        + numpy.multiply.outer(0.25*(1-xi[1]), quad[1, 0])
        + numpy.multiply.outer(0.25*(1+xi[1]), quad[1, 1])
        - numpy.multiply.outer(0.25*(1+xi[1]), quad[0, 1])
        ).T
    J1 = (
        - numpy.multiply.outer(0.25*(1-xi[0]), quad[0, 0])
        - numpy.multiply.outer(0.25*(1+xi[0]), quad[1, 0])
        + numpy.multiply.outer(0.25*(1+xi[0]), quad[1, 1])
        + numpy.multiply.outer(0.25*(1-xi[0]), quad[0, 1])
        ).T
    return J0[0]*J1[1] - J1[0]*J0[1]


def integrate(f, quad, scheme, sumfun=helpers.kahan_sum):
    x = transform(scheme.points.T, quad).T
    det = _get_detJ(scheme.points.T, quad)
    return sumfun(scheme.weights * f(x) * abs(det), axis=-1)
