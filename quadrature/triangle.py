# -*- coding: utf-8 -*-
#
import numpy


def volume(triangle):
    # It doesn't matter much which cross product we take for computing the
    # triangle volumes; deliberately take
    #
    #   <e0 x e1, e0 x e1> = <e0, e0> <e1, e1> - <e0, e1>^2.
    #
    e0 = triangle[1] - triangle[0]
    e1 = triangle[2] - triangle[1]
    e0_dot_e0 = numpy.dot(e0, e0)
    e0_dot_e1 = numpy.dot(e0, e1)
    e1_dot_e1 = numpy.dot(e1, e1)
    return 0.5 * numpy.sqrt(e0_dot_e0 * e1_dot_e1 - e0_dot_e1**2)


def _transform_to_unit_triangle(f, triangle):
    '''Transformation

      x = x0 * N0(xi, eta) + x1 * N1(xi, eta) + x2 * N2(xi, eta)

    with

      N0(xi, eta) = 1 - xi - eta,
      N1(xi, eta) = xi,
      N2(xi, eta) = eta.
    '''
    return lambda xi: f(
        + triangle[0] * (1.0 - xi[0] - xi[1])
        + triangle[1] * xi[0]
        + triangle[2] * xi[1]
        )


def integrate(f, triangle, rule):
    g = _transform_to_unit_triangle(f, triangle)
    out = 0.0
    for point, weight in zip(rule.points, rule.weights):
        out += weight * g(point)
    return volume(triangle) * out


class Centroid(object):
    def __init__(self):
        self.weights = [1.0]
        self.points = numpy.array([
            [1.0/3.0, 1.0/3.0]
            ])
        self.order = 1
        return


class Vertex(object):
    def __init__(self):
        self.weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
        self.points = numpy.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            ])
        self.order = 1
        return
