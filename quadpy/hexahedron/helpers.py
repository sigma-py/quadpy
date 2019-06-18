# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from ..ncube import transform, NCubeScheme, ncube_points as cube_points


class HexahedronScheme(NCubeScheme):

    def __init__(self, name, weights, points, degree):
        self.name = name
        self.weights = weights
        self.points = points
        self.degree = degree
        return

    def show(self, hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]), backend="vtk"):
        """Shows the quadrature points on a given hexahedron. The size of the balls
        around the points coincides with their weights."""
        edges = numpy.array(
            [
                [hexa[0, 0, 0], hexa[1, 0, 0]],
                [hexa[1, 0, 0], hexa[1, 1, 0]],
                [hexa[1, 1, 0], hexa[0, 1, 0]],
                [hexa[0, 1, 0], hexa[0, 0, 0]],
                #
                [hexa[0, 0, 1], hexa[1, 0, 1]],
                [hexa[1, 0, 1], hexa[1, 1, 1]],
                [hexa[1, 1, 1], hexa[0, 1, 1]],
                [hexa[0, 1, 1], hexa[0, 0, 1]],
                #
                [hexa[0, 0, 0], hexa[0, 0, 1]],
                [hexa[1, 0, 0], hexa[1, 0, 1]],
                [hexa[1, 1, 0], hexa[1, 1, 1]],
                [hexa[0, 1, 0], hexa[0, 1, 1]],
            ]
        )
        edges = numpy.moveaxis(edges, 1, 2)

        helpers.backend_to_function[backend](
            transform(self.points.T, hexa),
            self.weights,
            self.integrate(lambda x: 1.0, hexa),
            edges,
        )
        return


def z():
    return numpy.array([[0, 0, 0]])


def fs_r00(a):
    return numpy.array(
        [[+a, 0, 0], [0, +a, 0], [0, 0, +a], [-a, 0, 0], [0, -a, 0], [0, 0, -a]]
    )


def fs_rr0(a):
    return numpy.array(
        [
            [+a, +a, 0],
            [+a, 0, +a],
            [0, +a, +a],
            [+a, -a, 0],
            [+a, 0, -a],
            [0, +a, -a],
            [-a, +a, 0],
            [-a, 0, +a],
            [0, -a, +a],
            [-a, -a, 0],
            [-a, 0, -a],
            [0, -a, -a],
        ]
    )


def fs_rrs(a, b):
    return numpy.array(
        [
            [+a, +a, +b],
            [+a, +b, +a],
            [+b, +a, +a],
            [+a, -a, +b],
            [+a, +b, -a],
            [+b, +a, -a],
            [-a, +a, +b],
            [-a, +b, +a],
            [+b, -a, +a],
            [-a, -a, +b],
            [-a, +b, -a],
            [+b, -a, -a],
            [+a, +a, -b],
            [+a, -b, +a],
            [-b, +a, +a],
            [+a, -a, -b],
            [+a, -b, -a],
            [-b, +a, -a],
            [-a, +a, -b],
            [-a, -b, +a],
            [-b, -a, +a],
            [-a, -a, -b],
            [-a, -b, -a],
            [-b, -a, -a],
        ]
    )


def rss_pm(r, s):
    return numpy.array(
        [
            [+r, +s, +s],
            [+s, +r, +s],
            [+s, +s, +r],
            [-r, -s, -s],
            [-s, -r, -s],
            [-s, -s, -r],
        ]
    )


def pm_rrr(a):
    return numpy.array(
        [
            [+a, +a, +a],
            [-a, +a, +a],
            [+a, -a, +a],
            [-a, -a, +a],
            [+a, +a, -a],
            [-a, +a, -a],
            [+a, -a, -a],
            [-a, -a, -a],
        ]
    )
