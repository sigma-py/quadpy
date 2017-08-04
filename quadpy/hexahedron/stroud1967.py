# -*- coding: utf-8 -*-
#
import numpy

from .helpers import rss_pm, z

from ..helpers import untangle


class Stroud1967(object):
    '''
    A.H. Stroud,
    Some fifth degree integration formulas for symmetric regions II,
    Numerische Mathematik, Volume 9 Issue 5, April 1967, Pages 460-468
    <https://dx.doi.org/10.1007/BF02162160>.

    Analytic expression for all quantities are given in

    J.W. Peterson,
    Analytical Formulae for Two of A. H. Stroud's Quadrature Rules,
    Sep. 2009,
    <https://arxiv.org/pdf/0909.5106.pdf>.
    '''
    def __init__(self):
        self.degree = 5

        sqrt19 = numpy.sqrt(19.0)
        t = numpy.sqrt(71440.0 + 6802.0 * sqrt19)

        lmbd = numpy.sqrt((1919.0 - 148.0*sqrt19 + 4*t) / 3285.0)
        gmma = numpy.sqrt((1919.0 - 148.0*sqrt19 - 4*t) / 3285.0)
        xi = -numpy.sqrt((1121.0 + 74.0*sqrt19 - 2*t) / 3285.0)
        mu = +numpy.sqrt((1121.0 + 74.0*sqrt19 + 2*t) / 3285.0)

        B = 133225.0 / (260072.0 - 1520*sqrt19 + (133.0-37.0*sqrt19)*t)
        C = 133225.0 / (260072.0 - 1520*sqrt19 - (133.0-37.0*sqrt19)*t)

        data = [
            (32.0/19.0, z()),
            (B, rss_pm(lmbd, xi)),
            (C, rss_pm(gmma, mu)),
            ]

        self.points, self.weights = untangle(data)
        return
