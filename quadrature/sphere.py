# -*- coding: utf-8 -*-
#
import math
import numpy


def area(radius):
    return 4*numpy.pi*radius**2


# def show(sphere, scheme, circle_scale=1.0):
#     '''Shows the quadrature points on a given sphere. The size of the circles
#     around the points coincides with their weights.
#     '''
#     from matplotlib import pyplot as plt
#
#     plt.plot(sphere[:, 0], sphere[:, 1], '-k')
#     plt.plot(
#         [sphere[-1, 0], sphere[0, 0]],
#         [sphere[-1, 1], sphere[0, 1]],
#         '-k')
#
#     transformed_pts = \
#         + numpy.outer(
#             (1.0 - scheme.points[:, 0] - scheme.points[:, 1]),
#             sphere[0]
#             ) \
#         + numpy.outer(scheme.points[:, 0], sphere[1]) \
#         + numpy.outer(scheme.points[:, 1], sphere[2])
#
#     # plt.plot(transformed_pts[:, 0], transformed_pts[:, 1], 'or')
#     sphere_vol = volume(sphere)
#     for tp, weight in zip(transformed_pts, scheme.weights):
#         color = 'b' if weight >= 0 else 'r'
#         # highlight circle center
#         plt.plot([tp[0]], [tp[1]], '.' + color)
#         # plot circle
#         # scale the circle volume according to the weight
#         radius = circle_scale \
#             * numpy.sqrt(sphere_vol * abs(weight) / numpy.pi)
#         circ = plt.Circle((tp[0], tp[1]), radius, color=color, alpha=0.5)
#         plt.gcf().gca().add_artist(circ)
#
#     plt.axis('equal')
#     return


def integrate(f, midpoint, radius, rule):
    # w * f(x(xi)) * |det(J)|
    out = math.fsum([
        weight * radius**3 * f(radius*xi + midpoint)
        for xi, weight in zip(rule.points, rule.weights)
        ])
    return area(radius) * out


class Lebedev(object):
    '''
    Sphere integration schemes from

    Lebedev, V. I. (1976),
    Quadratures on a sphere,
    Zh. Vȳchisl. Mat. Mat. Fiz. 16 (2): 293–306.
    doi:10.1016/0041-5553(76)90100-2.

    <https://en.wikipedia.org/wiki/Lebedev_quadrature>
    <https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html>
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.concatenate([
                1.0/6.0 * numpy.ones(6),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                ])
            self.degree = 3
        elif index == 2:
            self.weights = numpy.concatenate([
                1.0/15.0 * numpy.ones(6),
                0.075 * numpy.ones(8),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a3(),
                ])
            self.degree = 5
        elif index == 3:
            self.weights = numpy.concatenate([
                1.0/21.0 * numpy.ones(6),
                4.0/105.0 * numpy.ones(12),
                9.0/280.0 * numpy.ones(8),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a2(),
                self.a3(),
                ])
            self.degree = 7
        else:
            raise ValueError('Illegal Lebedev index')

        return

    def a1(self):
        return numpy.array([
            [+1.0, 0.0, 0.0],
            [0.0, +1.0, 0.0],
            [0.0, 0.0, +1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            ])

    def a2(self):
        return numpy.array([
            [+1.0, +1.0, 0.0],
            [+1.0, 0.0, +1.0],
            [0.0, +1.0, +1.0],
            [-1.0, +1.0, 0.0],
            [-1.0, 0.0, +1.0],
            [0.0, -1.0, +1.0],
            [+1.0, -1.0, 0.0],
            [+1.0, 0.0, -1.0],
            [0.0, +1.0, -1.0],
            [-1.0, -1.0, 0.0],
            [-1.0, 0.0, -1.0],
            [0.0, -1.0, -1.0],
            ]) / numpy.sqrt(2.0)

    def a3(self):
        return numpy.array([
            [+1.0, +1.0, +1.0],
            [+1.0, +1.0, -1.0],
            [+1.0, -1.0, +1.0],
            [+1.0, -1.0, -1.0],
            [-1.0, +1.0, +1.0],
            [-1.0, +1.0, -1.0],
            [-1.0, -1.0, +1.0],
            [-1.0, -1.0, -1.0],
            ]) / numpy.sqrt(3.0)
