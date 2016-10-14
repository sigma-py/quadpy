# -*- coding: utf-8 -*-
#
import math
import numpy
from . import tetrahedron


def integrate(f, wedge, scheme):
    xi = scheme.points.T
    x = \
        + numpy.outer(wedge[0], 0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2])) \
        + numpy.outer(wedge[1], 0.5 * xi[0] * (1.0-xi[2])) \
        + numpy.outer(wedge[2], 0.5 * xi[1] * (1.0-xi[2])) \
        + numpy.outer(wedge[3], 0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2])) \
        + numpy.outer(wedge[4], 0.5 * xi[0] * (1.0+xi[2])) \
        + numpy.outer(wedge[5], 0.5 * xi[1] * (1.0+xi[2]))
    J0 = \
        - numpy.outer(wedge[0], 0.5*(1.0 - xi[2])) \
        + numpy.outer(wedge[2], 0.5*(1.0 - xi[2])) \
        - numpy.outer(wedge[3], 0.5*(1.0 + xi[2])) \
        + numpy.outer(wedge[5], 0.5*(1.0 + xi[2]))
    J1 = \
        - numpy.outer(wedge[0], 0.5*(1.0 - xi[2])) \
        + numpy.outer(wedge[1], 0.5*(1.0 - xi[2])) \
        - numpy.outer(wedge[3], 0.5*(1.0 + xi[2])) \
        + numpy.outer(wedge[4], 0.5*(1.0 + xi[2]))
    J2 = \
        - numpy.outer(wedge[0], 0.5 * (1.0-xi[0]-xi[1])) \
        - numpy.outer(wedge[1], 0.5 * xi[0]) \
        - numpy.outer(wedge[2], 0.5 * xi[1]) \
        + numpy.outer(wedge[3], 0.5 * (1.0-xi[0]-xi[1])) \
        + numpy.outer(wedge[4], 0.5 * xi[0]) \
        + numpy.outer(wedge[5], 0.5 * xi[1])
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]

    return math.fsum(scheme.weights * f(x).T * abs(det))


class Felippa(object):
    '''
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.

    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_wedge/quadrature_rules_wedge.html>
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.array([1.0])
            self.points = numpy.array([
                [1.0/3.0, 1.0/3.0, 0.0]
                ])
            self.degree = 1
        elif index == 2:
            self.weights = 1.0/6.0 * numpy.ones(6)
            self.points = numpy.array([
                [2.0/3.0, 1.0/6.0, -numpy.sqrt(1.0/3.0)],
                [1.0/6.0, 2.0/3.0, -numpy.sqrt(1.0/3.0)],
                [1.0/6.0, 1.0/6.0, -numpy.sqrt(1.0/3.0)],
                [2.0/3.0, 1.0/6.0, +numpy.sqrt(1.0/3.0)],
                [1.0/6.0, 2.0/3.0, +numpy.sqrt(1.0/3.0)],
                [1.0/6.0, 1.0/6.0, +numpy.sqrt(1.0/3.0)],
                ])
            self.degree = 2
        else:
            raise ValueError('Illegal Felippa order')

        return
