# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle


class Felippa(object):
    """
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.

    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_wedge/quadrature_rules_wedge.html>
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        if index == 1:
            self.degree = 1
            data = [(1, _s3(symbolic))]
        elif index == 2:
            self.degree = 2
            data = [(frac(1, 6), _s21_z(frac(1, 6), sqrt(frac(1, 3))))]
        elif index == 3:
            self.degree = 2
            data = [(frac(1, 6), _s21_z(frac(1, 2), sqrt(frac(1, 3))))]
        elif index == 4:
            self.degree = 4
            # roots of  135 x^4 - 240 x^3 + 120 x^2 - 20 x + 1
            a1, a2 = [
                (40 - 5 * sqrt(10) - i * sqrt(950 - 220 * sqrt(10))) / 90
                for i in [+1, -1]
            ]
            data = [
                (0.6205044157722541E-01, _s21_z(a2, sqrt(frac(3, 5)))),
                (0.3054215101536719E-01, _s21_z(a1, sqrt(frac(3, 5)))),
                (0.9928070652356065E-01, _s21(a2)),
                (0.4886744162458750E-01, _s21(a1)),
            ]
        elif index == 5:
            self.degree = 5
            a1, a2 = [(6 - i * sqrt(15)) / 21 for i in [+1, -1]]
            data = [
                (0.3498310570689643E-01, _s21_z(a1, sqrt(frac(3, 5)))),
                (0.3677615355236283E-01, _s21_z(a2, sqrt(frac(3, 5)))),
                (frac(1, 16), _s3_z(sqrt(frac(3, 5)), symbolic)),
                (0.5597296913103428E-01, _s21(a1)),
                (0.5884184568378053E-01, _s21(a2)),
                (frac(1, 10), _s3(symbolic)),
            ]
        else:
            assert index == 6
            self.degree = 6

            data = [
                (
                    0.8843323515718317E-02,
                    _s21_z(0.6308901449150223E-01, -0.8611363115940526),
                ),
                (
                    0.2031233592848984E-01,
                    _s21_z(0.2492867451709104, -0.8611363115940526),
                ),
                (
                    0.1441007403935041E-01,
                    _s111_z(
                        0.5314504984481695E-01, 0.3103524510337844, 0.8611363115940526
                    ),
                ),
                (
                    0.1657912966938509E-01,
                    _s21_z(0.6308901449150223E-01, 0.3399810435848563),
                ),
                (
                    0.3808080193469984E-01,
                    _s21_z(0.2492867451709104, 0.3399810435848563),
                ),
                (
                    0.2701546376983638E-01,
                    _s111_z(
                        0.5314504984481695E-01, 0.3103524510337844, 0.3399810435848563
                    ),
                ),
            ]

        self.points, self.weights = untangle(data)
        return


def _s3(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return [[frac(1, 3), frac(1, 3), 0]]


def _s3_z(z, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return [[frac(1, 3), frac(1, 3), +z], [frac(1, 3), frac(1, 3), -z]]


def _s21(a):
    b = 1 - 2 * a
    return [[a, b, 0], [b, a, 0], [a, a, 0]]


def _s21_z(a, z):
    b = 1 - 2 * a
    return [[a, b, +z], [b, a, +z], [a, a, +z], [a, b, -z], [b, a, -z], [a, a, -z]]


def _s111_z(a, b, z):
    c = 1 - a - b
    return [
        [b, c, +z],
        [a, b, +z],
        [c, a, +z],
        [c, b, +z],
        [a, c, +z],
        [b, a, +z],
        [b, c, -z],
        [a, b, -z],
        [c, a, -z],
        [c, b, -z],
        [a, c, -z],
        [b, a, -z],
    ]
