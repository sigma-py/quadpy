# -*- coding: utf-8 -*-
#
from sympy import Rational as fr, sqrt

from .helpers import rss_pm, z
from ..helpers import untangle


class Stroud1967(object):
    '''
    A.H. Stroud,
    Some fifth degree integration formulas for symmetric regions II,
    Numerische Mathematik, Volume 9 Issue 5, April 1967, Pages 460-468
    <https://doi.org/10.1007/BF02162160>.

    Analytic expression for all quantities are given in

    J.W. Peterson,
    Analytical Formulae for Two of A. H. Stroud's Quadrature Rules,
    Sep. 2009,
    <https://arxiv.org/pdf/0909.5106.pdf>.
    '''
    def __init__(self):
        self.degree = 5

        sqrt19 = sqrt(19)
        t = sqrt(71440 + 6802 * sqrt19)

        lmbd, gmma = [
            sqrt((1919 - 148*sqrt19 + i * 4*t) / 3285)
            for i in [+1, -1]
            ]
        xi, mu = [
            -sqrt((1121 + 74*sqrt19 - i * 2*t) / 3285)
            for i in [+1, -1]
            ]
        mu *= -1

        B, C = [
            133225 / (260072 - 1520*sqrt19 + i*(133 - 37*sqrt19)*t)
            for i in [+1, -1]
            ]

        data = [
            (fr(32, 19), z()),
            (B, rss_pm(lmbd, xi)),
            (C, rss_pm(gmma, mu)),
            ]

        self.points, self.weights = untangle(data)
        return
