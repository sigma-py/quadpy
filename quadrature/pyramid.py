# -*- coding: utf-8 -*-
#
import math
import numpy


def integrate(f, pyra, scheme):
    xi = scheme.points.T
    x = \
        + numpy.outer(pyra[0], 0.125*(1.0-xi[0])*(1.0-xi[1])*(1-xi[2])) \
        + numpy.outer(pyra[1], 0.125*(1.0+xi[0])*(1.0-xi[1])*(1-xi[2])) \
        + numpy.outer(pyra[2], 0.125*(1.0+xi[0])*(1.0+xi[1])*(1-xi[2])) \
        + numpy.outer(pyra[3], 0.125*(1.0-xi[0])*(1.0+xi[1])*(1-xi[2])) \
        + numpy.outer(pyra[4], 0.50*(1.0+xi[2]))
    J0 = \
        - numpy.outer(pyra[0], 0.125*(1.0-xi[1])*(1-xi[2])) \
        + numpy.outer(pyra[1], 0.125*(1.0-xi[1])*(1-xi[2])) \
        + numpy.outer(pyra[2], 0.125*(1.0+xi[1])*(1-xi[2])) \
        - numpy.outer(pyra[3], 0.125*(1.0+xi[1])*(1-xi[2]))
    J1 = \
        - numpy.outer(pyra[0], 0.125*(1.0-xi[0])*(1-xi[2])) \
        - numpy.outer(pyra[1], 0.125*(1.0+xi[0])*(1-xi[2])) \
        + numpy.outer(pyra[2], 0.125*(1.0+xi[0])*(1-xi[2])) \
        + numpy.outer(pyra[3], 0.125*(1.0-xi[0])*(1-xi[2]))
    J2 = \
        - numpy.outer(pyra[0], 0.125*(1.0-xi[0])*(1.0-xi[1])) \
        - numpy.outer(pyra[1], 0.125*(1.0+xi[0])*(1.0-xi[1])) \
        - numpy.outer(pyra[2], 0.125*(1.0+xi[0])*(1.0+xi[1])) \
        - numpy.outer(pyra[3], 0.125*(1.0-xi[0])*(1.0+xi[1])) \
        + numpy.outer(pyra[4], 0.50*numpy.ones(1))
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]

    return math.fsum(scheme.weights * f(x).T * abs(det))


class Felippa(object):
    '''
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.

    Abstract:
    This paper presents a set of Mathematica modules that organizes numerical
    integration rules considered useful for finite element work. Seven regions
    are considered: line segments, triangles, quadrilaterals, tetrahedral,
    wedges, pyramids and hexahedra. Information can be returned in
    floating-point (numerical) form, or in exact symbolic form. The latter is
    useful for computer-algebra aided FEM work that carries along symbolic
    variables. A few quadrature rules were extracted from sources in the FEM
    and computational mathematics literature, and placed in symbolic form using
    Mathematica to generate own code. A larger class of formulas, previously
    known only numerically, were directly obtained through symbolic
    computations. Some unpublished non-product rules for pyramid regions were
    found and included in the collection. For certain regions: quadrilaterals,
    wedges and hexahedra, only product rules were included to economize
    programming. The collection embodies most FEM-useful formulas of low and
    moderate order for the seven regions noted above. Some gaps as regard
    region geometries and omission of non-product rules are noted in the
    conclusions. The collection may be used “as is” in support of symbolic FEM
    work thus avoiding contamination with floating arithmetic that precludes
    simplification. It can also be used as generator for low-level
    floating-point code modules in Fortran or C. Floating point accuracy can be
    selected arbitrarily. No similar modular collection applicable to a range
    of FEM work, whether symbolic or numeric, has been published before.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.array([128.0/27.0])
            self.points = numpy.array([
                [0.0, 0.0, -0.5],
                ])
            self.degree = 1
        elif index == 2:
            self.weights = numpy.concatenate([
                0.81 * numpy.ones(4),
                125.0/27.0 * numpy.ones(1)
                ])
            self.points = numpy.concatenate([
                self._s4(8 * numpy.sqrt(2.0/15.0) / 5, -2.0/3.0),
                numpy.array([[0.0, 0.0, 0.4]]),
                ])
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                504.0/625.0 * numpy.ones(4),
                576.0/625.0 * numpy.ones(1),
                64.0/15.0 * numpy.ones(1),
                ])
            self.points = numpy.concatenate([
                self._s4(numpy.sqrt(12.0/35.0), -2.0/3.0),
                numpy.array([[0.0, 0.0, 1.0/6.0]]),
                numpy.array([[0.0, 0.0, 0.5]]),
                ])
            self.degree = 2
        elif index == 4:
            w1 = 5 * (68.0 + 5*numpy.sqrt(10.0)) / 432.0
            w2 = 85.0/54.0 - w1
            self.weights = numpy.concatenate([
                w1 * numpy.ones(4),
                w2 * numpy.ones(4),
                ])
            g1 = numpy.sqrt(1.0/3.0)
            g2 = (2*numpy.sqrt(10)-5) / 15.0
            self.points = numpy.concatenate([
                self._s4(g1, g2),
                self._s4(g1, -2.0/3.0 - g2),
                ])
            self.degree = 3
        elif index == 5:
            w1 = (11764.0 - 461.0*numpy.sqrt(51.0)) / 15300.0
            w2 = 346.0 / 225.0 - w1
            self.weights = numpy.concatenate([
                w1 * numpy.ones(4),
                w2 * numpy.ones(4),
                ])
            g1 = numpy.sqrt(2.0/15.0 * (573 - 2*numpy.sqrt(51))) / 15.0
            g2 = numpy.sqrt(2.0/15.0 * (573 + 2*numpy.sqrt(51))) / 15.0
            g3 = -(2*numpy.sqrt(51.0) + 13) / 35.0
            g4 = +(2*numpy.sqrt(51.0) - 13) / 35.0
            self.points = numpy.concatenate([
                self._s4(g1, g3),
                self._s4(g2, g4),
                ])
            self.degree = 2
        elif index == 6:
            w1 = 7.0*(11472415.0 - 70057.0*numpy.sqrt(2865.0)) / 130739500.0
            w2 = 84091.0/68450.0 - w1
            self.weights = numpy.concatenate([
                w1 * numpy.ones(4),
                w2 * numpy.ones(4),
                3.6 * numpy.ones(1),
                ])

            g1 = 8 * numpy.sqrt(
                (573 + 5*numpy.sqrt(2865.0))
                / (109825 + 969*numpy.sqrt(2865.0))
                )
            g2 = numpy.sqrt(2*(8025 + numpy.sqrt(2865.0)) / 35.0) / 37.0
            g3 = -(+87 + numpy.sqrt(2865.0)) / 168.0
            g4 = +(-87 + numpy.sqrt(2865.0)) / 168.0
            self.points = numpy.concatenate([
                self._s4(g1, g3),
                self._s4(g2, g4),
                numpy.array([[0.0, 0.0, 2.0/3.0]])
                ])
            self.degree = 2
        else:
            raise ValueError('Illegal Felippa order')

        return

    def _s4(self, a, z):
        return numpy.array([
            [+a, +a, z],
            [-a, +a, z],
            [+a, -a, z],
            [-a, -a, z],
            ])
