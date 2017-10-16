# -*- coding: utf-8 -*-
#
import numpy
from sympy import Rational as fr, sqrt

from .helpers import _s4, _s4_0
from ..helpers import untangle


class Felippa(object):
    '''
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890,
    <http://doi.org/10.1108/02644400410554362>.

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
        wg9 = numpy.array([
            fr(64, 81),
            fr(40, 81),
            fr(25, 81),
            ])

        if index == 1:
            self.degree = 1
            data = [
                (fr(128, 27), numpy.array([[0, 0, -fr(1, 2)]])),
                ]
        elif index == 2:
            self.degree = 2
            data = [
                (fr(81, 100), _s4(8 * sqrt(fr(2, 15)) / 5, -fr(2, 3))),
                (fr(125, 27), numpy.array([[0, 0, fr(2, 5)]])),
                ]
        elif index == 3:
            self.degree = 2
            data = [
                (fr(504, 625), _s4(sqrt(fr(12, 35)), -fr(2, 3))),
                (fr(576, 625), numpy.array([[0, 0, fr(1, 6)]])),
                (fr(64, 15), numpy.array([[0, 0, fr(1, 2)]])),
                ]
        elif index == 4:
            self.degree = 3
            w1 = 5 * (68 + 5*sqrt(10)) / 432
            w2 = fr(85, 54) - w1
            g1 = sqrt(fr(1, 3))
            g2 = (2*sqrt(10)-5) / 15
            data = [
                (w1, _s4(g1, g2)),
                (w2, _s4(g1, -fr(2, 3) - g2)),
                ]
        elif index == 5:
            self.degree = 2
            w1 = (11764 - 461*sqrt(51)) / 15300
            w2 = fr(346, 225) - w1
            g1, g2 = [
                sqrt(fr(2, 15) * (573 - i * 2*sqrt(51))) / 15
                for i in [+1, -1]
                ]
            g3, g4 = [
                -i * (2*sqrt(51) + i*13) / 35
                for i in [+1, -1]
                ]
            data = [
                (w1, _s4(g1, g3)),
                (w2, _s4(g2, g4)),
                ]
        elif index == 6:
            self.degree = 2
            w1 = 7*(11472415 - 70057*sqrt(2865)) / 130739500
            w2 = fr(84091, 68450) - w1

            g1 = 8 * sqrt(
                (573 + 5*sqrt(2865))
                / (109825 + 969*sqrt(2865))
                )
            g2 = sqrt(2*(8025 + sqrt(2865)) / 35) / 37
            g3, g4 = [
                -i * (+i*87 + sqrt(2865)) / 168
                for i in [+1, -1]
                ]

            data = [
                (w1, _s4(g1, g3)),
                (w2, _s4(g2, g4)),
                (fr(18, 5), numpy.array([[0, 0, fr(2, 3)]])),
                ]
        elif index == 7:
            self.degree = 2
            w1 = fr(170569, 331200)
            w2 = fr(276710106577408, 1075923777052725)
            w3 = fr(12827693806929, 30577384040000)
            w4 = fr(10663383340655070643544192, 4310170528879365193704375)
            g1 = 7 * sqrt(fr(35, 59)) / 8
            g2 = 224 * sqrt(fr(336633710, 33088740423)) / 37
            g3 = sqrt(fr(37043, 35)) / 56
            g4 = -fr(127, 153)
            g5 = fr(1490761, 2842826)
            data = [
                (w1, _s4(g1, -fr(1, 7))),
                (w2, _s4_0(g2, -fr(9, 28))),
                (w3, _s4(g3, g4)),
                (w4, numpy.array([[0, 0, g5]])),
                ]
        elif index == 8:
            self.degree = 3
            w1 = 5 * (68 + 5*sqrt(10)) / 432
            w2 = fr(85, 54) - w1
            g1 = sqrt(fr(3, 5))
            g2 = 1 - 2*(10 - sqrt(10)) / 15
            g3 = -fr(2, 3) - g2
            data = [
                (w1*wg9[2], _s4(g1, g2)),
                (w1*wg9[1], _s4_0(g1, g2)),
                (w1*wg9[0], numpy.array([[0, 0, g2]])),
                (w2*wg9[2], _s4(g1, g3)),
                (w2*wg9[1], _s4_0(g1, g3)),
                (w2*wg9[0], numpy.array([[0, 0, g3]])),
                ]
        else:
            assert index == 9
            self.degree = 5
            g1 = sqrt(fr(3, 5))
            g3 = -0.854011951853700535688324041975993416
            g4 = -0.305992467923296230556472913192103090
            g5 = +0.410004419776996766244796955168096505
            w1 = fr(4, 15)*(4+5*(g4+g5)+10*g4*g5)/((g3-g4)*(g3-g5)*(1-g3)**2)
            w2 = fr(4, 15)*(4+5*(g3+g5)+10*g3*g5)/((g3-g4)*(g5-g4)*(1-g4)**2)
            w3 = fr(4, 15)*(4+5*(g3+g4)+10*g3*g4)/((g3-g5)*(g4-g5)*(1-g5)**2)
            data = [
                (w1*wg9[2], _s4(g1, g3)),
                (w1*wg9[1], _s4_0(g1, g3)),
                (w1*wg9[0], numpy.array([[0, 0, g3]])),
                (w2*wg9[2], _s4(g1, g4)),
                (w2*wg9[1], _s4_0(g1, g4)),
                (w2*wg9[0], numpy.array([[0, 0, g4]])),
                (w3*wg9[2], _s4(g1, g5)),
                (w3*wg9[1], _s4_0(g1, g5)),
                (w3*wg9[0], numpy.array([[0, 0, g5]])),
                ]

        self.points, self.weights = untangle(data)
        return
