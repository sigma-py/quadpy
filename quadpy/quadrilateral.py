# -*- coding: utf-8 -*-
#
import numpy

from . import line_segment
from . import helpers


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


class Product(object):
    def __init__(self, scheme1d):
        self.schemes = \
            scheme1d if isinstance(scheme1d, list) \
            else 2 * [scheme1d]

        self.weights = numpy.outer(
            self.schemes[0].weights, self.schemes[1].weights
            ).flatten()
        self.points = numpy.dstack(numpy.meshgrid(
            self.schemes[0].points, self.schemes[1].points
            )).reshape(-1, 2)
        self.degree = min([s.degree for s in self.schemes])
        return


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.

    <https://people.sc.fsu.edu/~jburkardt/m_src/stroud/square_unit_set.m>
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = [4.0]
            self.points = [[0.0, 0.0]]
            self.degree = 1
        elif index == 2:
            self.weights = 4 * [1.0]
            self.points = _symm_s(1.0/numpy.sqrt(3.0))
            self.degree = 3
        elif index == 3:
            self.weights = (
                [64.0/81.0] +
                4 * [25.0/81.0] +
                4 * [40.0/81.0]
                )
            self.points = (
                [[0.0, 0.0]] +
                _symm_s(numpy.sqrt(0.6)) +
                _symm_r_0(numpy.sqrt(0.6))
                )
            self.degree = 5
        elif index == 4:
            # Stroud number C2:7-1.
            r = numpy.sqrt(6.0 / 7.0)
            c = 3.0 * numpy.sqrt(583.0)
            s = numpy.sqrt((114.0 - c) / 287.0)
            t = numpy.sqrt((114.0 + c) / 287.0)
            w1 = 4.0 * 49.0 / 810.0
            w2 = 4.0 * (178981.0 + 923.0 * c) / 1888920.0
            w3 = 4.0 * (178981.0 - 923.0 * c) / 1888920.0
            #
            self.weights = (
                4 * [w1] +
                4 * [w2] +
                4 * [w3]
                )
            self.points = (
                _symm_r_0(r) +
                _symm_s(s) +
                _symm_s(t)
                )
            self.degree = 7
        elif index == 5:
            # Stroud number C2:7-3.
            r = numpy.sqrt(12.0 / 35.0)
            c = 3.0 * numpy.sqrt(186.0)
            s = numpy.sqrt((93.0 + c) / 155.0)
            t = numpy.sqrt((93.0 - c) / 155.0)
            w1 = 8.0 / 162.0
            w2 = 98.0 / 162.0
            w3 = 31.0 / 162.0
            self.weights = (
                [w1] +
                4 * [w2] +
                8 * [w3]
                )
            self.points = (
                [[0.0, 0.0]] +
                _symm_r_0(r) +
                _symm_s_t(s, t)
                )
            self.degree = 7
        else:
            assert index == 6
            scheme1d = line_segment.GaussLegendre(8)
            self.weights = numpy.outer(
                scheme1d.weights, scheme1d.weights
                ).flatten()
            self.points = numpy.dstack(numpy.meshgrid(
                scheme1d.points, scheme1d.points
                )).reshape(-1, 2)
            assert len(self.points) == 64
            self.degree = 15

        self.weights = numpy.array(self.weights)
        self.points = numpy.array(self.points)
        return


def _symm_r_0(r):
    return [
        [+r, 0.0],
        [-r, 0.0],
        [0.0, +r],
        [0.0, -r],
        ]


def _symm_s(s):
    return [
        [+s, +s],
        [-s, +s],
        [+s, -s],
        [-s, -s],
        ]


def _symm_s_t(s, t):
    return [
        [+s, +t],
        [-s, +t],
        [+s, -t],
        [-s, -t],
        [+t, +s],
        [-t, +s],
        [+t, -s],
        [-t, -s],
        ]
