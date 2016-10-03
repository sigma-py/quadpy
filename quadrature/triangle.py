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
        self.degree = 1
        return


class Vertex(object):
    def __init__(self):
        self.weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
        self.points = numpy.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            ])
        self.degree = 1
        return


class SevenPoint(object):
    def __init__(self):
        self.weights = [
            0.45,
            0.05,
            0.05,
            0.05,
            2.0 / 15.0,
            2.0 / 15.0,
            2.0 / 15.0
            ]
        self.points = numpy.array([
            [1.0/3.0, 1.0/3.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.5],
            [0.5, 0.0],
            [0.5, 0.5],
            ])
        self.degree = 3
        return


class Gauss4x4(object):
    # https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    def __init__(self):
        self.weights = [
            0.04713673637581137,
            0.07077613579259895,
            0.04516809856187617,
            0.01084645180365496,
            0.08837017702418863,
            0.1326884322074010,
            0.08467944903812383,
            0.02033451909634504,
            0.08837017702418863,
            0.1326884322074010,
            0.08467944903812383,
            0.02033451909634504,
            0.04713673637581137,
            0.07077613579259895,
            0.04516809856187617,
            0.01084645180365496,
            ]

        self.points = numpy.array([
            [0.0571041961, 0.06546699455602246],
            [0.2768430136, 0.05021012321401679],
            [0.5835904324, 0.02891208422223085],
            [0.8602401357, 0.009703785123906346],
            [0.0571041961, 0.3111645522491480],
            [0.2768430136, 0.2386486597440242],
            [0.5835904324, 0.1374191041243166],
            [0.8602401357, 0.04612207989200404],
            [0.0571041961, 0.6317312516508520],
            [0.2768430136, 0.4845083266559759],
            [0.5835904324, 0.2789904634756834],
            [0.8602401357, 0.09363778440799593],
            [0.0571041961, 0.8774288093439775],
            [0.2768430136, 0.6729468631859832],
            [0.5835904324, 0.3874974833777692],
            [0.8602401357, 0.1300560791760936],
            ])

        self.degree = 7
        return
