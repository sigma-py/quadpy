# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
import numpy


class WilliamsShunnJameson(object):
    '''
    D.M. Williams, L. Shunn, A. Jameson,
    Symmetric quadrature rules for simplexes based on sphere close packed
    lattice arrangements,
    Journal of Computational and Applied Mathematics,
    266 (2014) 18–38.

    Abstract:
    Sphere close packed (SCP) lattice arrangements of points are well-suited
    for formulating symmetric quadrature rules on simplexes, as they are
    symmetric under affine transformations of the simplex unto itself in 2D and
    3D. As a result, SCP lattice arrangements have been utilized to formulate
    symmetric quadrature rules with Np = 1, 4, 10, 20, 35, and 56 points on the
    3-simplex (Shunn and Ham, 2012). In what follows, the work on the 3-simplex
    is extended, and SCP lattices are employed to identify symmetric quadrature
    rules with Np = 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, and 66 points on the
    2-simplex and Np = 84 points on the 3-simplex. These rules are found to be
    capable of exactly integrating polynomials of up to degree 17 in 2D and up
    to degree 9 in 3D.
    '''
    def __init__(self, index):
        self.name = 'WSJ(%d)' % index
        if index == 1:
            self.weights = [1.0]
            bary = _s3()
            self.degree = 1
        elif index == 2:
            self.weights = 3 * [1.0 / 3.0]
            bary = _s21(1.0 / 6.0)
            self.degree = 2
        elif index == 3:
            self.weights = (
                3 * [0.109951743655333] +
                3 * [0.223381589678000]
                )
            bary = numpy.concatenate([
                _s21(0.091576213509780),
                _s21(0.445948490915964)
                ])
            self.degree = 4
        else:
            assert index == 4
            self.weights = (
                3 * [0.041955512996649] +
                6 * [0.112098412070887] +
                [0.201542988584730]
                )
            bary = numpy.concatenate([
                _s21(0.055564052669793),
                _s111(0.295533711735893, 0.634210747745723),
                _s3(),
                ])
            self.degree = 5

        self.weights = numpy.array(self.weights)
        bary = numpy.array(bary)
        self.points = bary[:, [1, 2]]
        return
