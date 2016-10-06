# -*- coding: utf-8 -*-
#
import numpy


def volume(quad):
    # It doesn't matter much which cross product we take for computing the
    # quad volumes; deliberately take
    #
    #   <e0 x e1, e0 x e1> = <e0, e0> <e1, e1> - <e0, e1>^2.
    #
    e0 = quad[1] - quad[0]
    e1 = quad[2] - quad[1]
    e2 = quad[3] - quad[1]
    e0_dot_e0 = numpy.dot(e0, e0)
    e0_dot_e1 = numpy.dot(e0, e1)
    e1_dot_e1 = numpy.dot(e1, e1)
    e1_dot_e2 = numpy.dot(e1, e2)
    e2_dot_e2 = numpy.dot(e2, e2)
    return 0.5 * (
        numpy.sqrt(e0_dot_e0 * e1_dot_e1 - e0_dot_e1**2) +
        numpy.sqrt(e1_dot_e1 * e2_dot_e2 - e1_dot_e2**2)
        )


def show(quad, scheme, circle_scale=1.0):
    '''Shows the quadrature points on a given quad. The size of the circles
    around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt

    plt.plot(quad[:, 0], quad[:, 1], '-k')
    plt.plot(
        [quad[-1, 0], quad[0, 0]],
        [quad[-1, 1], quad[0, 1]],
        '-k')

    xi = scheme.points[:, 0]
    eta = scheme.points[:, 1]
    transformed_pts = \
        + numpy.outer(0.25 * (1.0 - xi)*(1.0 - eta), quad[0]) \
        + numpy.outer(0.25 * (1.0 + xi)*(1.0 - eta), quad[1]) \
        + numpy.outer(0.25 * (1.0 + xi)*(1.0 + eta), quad[2]) \
        + numpy.outer(0.25 * (1.0 - xi)*(1.0 + eta), quad[3])

    quad_vol = volume(quad)
    for tp, weight in zip(transformed_pts, scheme.weights):
        color = 'b' if weight >= 0 else 'r'
        # highlight circle center
        plt.plot([tp[0]], [tp[1]], '.' + color)
        # plot circle
        # scale the circle volume according to the weight
        radius = circle_scale * \
            numpy.sqrt(quad_vol * abs(weight)/sum(scheme.weights) / numpy.pi)
        circ = plt.Circle((tp[0], tp[1]), radius, color=color, alpha=0.5)
        plt.gcf().gca().add_artist(circ)

    plt.axis('equal')
    return


def integrate(f, quad, scheme):
    out = 0.0
    for xi, weight in zip(scheme.points, scheme.weights):
        x = \
            + quad[0] * 0.25*(1.0 - xi[0])*(1.0 - xi[1]) \
            + quad[1] * 0.25*(1.0 + xi[0])*(1.0 - xi[1]) \
            + quad[2] * 0.25*(1.0 + xi[0])*(1.0 + xi[1]) \
            + quad[3] * 0.25*(1.0 - xi[0])*(1.0 + xi[1])
        J0 = \
            - quad[0] * 0.25*(1-xi[1]) \
            + quad[1] * 0.25*(1-xi[1]) \
            + quad[2] * 0.25*(1+xi[1]) \
            - quad[3] * 0.25*(1+xi[1])
        J1 = \
            - quad[0] * 0.25*(1-xi[0]) \
            - quad[1] * 0.25*(1+xi[0]) \
            + quad[2] * 0.25*(1+xi[0]) \
            + quad[3] * 0.25*(1-xi[0])

        det_J = J0[0]*J1[1] - J1[0]*J0[1]
        out += weight * f(x) * abs(det_J)
    return out


class From1d(object):
    def __init__(self, scheme1d):
        self.weights = numpy.outer(
            scheme1d.weights, scheme1d.weights
            ).flatten()
        self.points = numpy.dstack(numpy.meshgrid(
            scheme1d.points, scheme1d.points
            )).reshape(-1, 2)
        self.degree = scheme1d.degree
        return
