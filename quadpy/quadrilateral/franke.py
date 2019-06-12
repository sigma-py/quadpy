# -*- coding: utf-8 -*-
#
"""
Richard Franke,
Obtaining cubatures for rectangles and other planar regions by using orthogonal
polynomials,
Math. Comp. 25 (1971), 803-817,
<https://doi.org/10.1090/S0025-5718-1971-0300440-5>.
"""
from __future__ import division

import numpy
import sympy

from .helpers import pmx, pmy, pm2, concat


class Franke1(object):
    def __init__(self, lmbda, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "Franke(1, {})".format(lmbda)

        assert -frac(9, 5) <= lmbda <= frac(9, 4)

        self.degree = 5

        a = sqrt(frac(9 + 5 * lmbda, 15))
        b = sqrt(frac(9 - 4 * lmbda, 15))
        c = sqrt(frac(3, 5))

        self.points, self.weights = concat(
            ([[0, 0]], [frac(16 * (4 + 5 * lmbda), 9 * (9 + 5 * lmbda))]),
            pmx([frac(40, 9 * (9 + 5 * lmbda)), a]),
            pm2([frac(25, 9 * (9 - 4 * lmbda)), b, c]),
            pmy([frac(40 * (1 - lmbda), 9 * (9 - 4 * lmbda)), c]),
        )
        return


class Franke2a(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 7

        a = sqrt(frac(15 + 2 * sqrt(30), 35))
        b = sqrt(frac(15 - 2 * sqrt(30), 35))

        self.points, self.weights = concat(
            pm2(
                [0.437841520872291e-1, 0.105784012371275e1, a],
                [0.362302863812526, 0.774596669241483, b],
                [0.304070693050225, 0.469253522127911, a],
            ),
            pmy([0.579684582100041, b]),
        )
        return


class Franke2b(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 7

        a = sqrt(frac(15 + 2 * sqrt(30), 35))
        b = sqrt(frac(15 - 2 * sqrt(30), 35))

        self.points, self.weights = concat(
            pm2(
                [0.193252691743030, 0.774596669241483, a],
                [0.169049921219002, 0.915060523380880, b],
                [0.483095233643544, 0.396191039748320, b],
            ),
            pmy([0.309204306788848, a]),
        )
        return


class Franke3a(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 9

        a = sqrt(frac(5, 9) + frac(2, 63) * sqrt(70))
        b = sqrt(frac(5, 9) - frac(2, 63) * sqrt(70))

        self.points, self.weights = concat(
            pm2(
                [0.705065140564012e-1, 0.845927799771709, a],
                [0.721121511007611e-1, 0.628901636732253, a],
                [0.971492736037507e-1, 0.959681421214621, b],
                [0.368549048677049, 0.436030596273468, b],
            ),
            pmx([0.316049382716049, 0.774596669241483]),
            pmy([0.188616439798053, a], [0.258606964371341e-1, b]),
            ([[0, 0]], [0.505679012345679]),
        )
        return


Franke = {"1": Franke1, "2a": Franke2a, "2b": Franke2b, "3a": Franke3a}
