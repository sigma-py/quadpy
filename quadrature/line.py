# -*- coding: utf-8 -*-
#
import numpy


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


def integrate(f, a, b, rule):
    out = 0.0
    for point, weight in zip(rule.points, rule.weights):
        out += weight * f(0.5 * (point + 1) * (b-a) + a)
    return (b - a) * out


class Midpoint(object):
    def __init__(self):
        self.weights = [1.0]
        self.points = numpy.array([
            0.0
            ])
        self.degree = 1
        return
