# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _s3
from ..helpers import untangle


class LiuVinokur(object):
    """
    Y. Liu and M. Vinokur,
    Exact Integrations of Polynomials and Symmetric Quadrature Formulas over
    Arbitrary Polyhedral Grids,
    Journal of Computational Physics, 140, 122–147 (1998).
    DOI: 10.1006/jcph.1998.5884,
    <https://doi.org/10.1006/jcph.1998.5884>.
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.frac if symbolic else lambda x, y: x / y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.name = "LV(%d)" % index
        if index == 1:
            self.degree = 1
            data = [(1, _s3(symbolic))]
        elif index == 2:
            self.degree = 1
            data = [(frac(1, 3), _r_alpha(1))]
        elif index == 3:
            self.degree = 2
            data = [(frac(1, 3), _r_alpha(-frac(1, 2)))]
        elif index == 4:
            self.degree = 2
            data = [(frac(3, 4), _s3(symbolic)), (frac(1, 12), _r_alpha(1))]
        elif index == 5:
            self.degree = 3
            data = [
                (-frac(9, 16), _s3(symbolic)),
                # ERR Wrongly specified in the article as 25 (instead of 2/5).
                (frac(25, 48), _r_alpha(frac(2, 5))),
            ]
        elif index == 6:
            self.degree = 3
            sqrt21 = sqrt(21)
            data = [
                ((1 + sqrt21) / 120, _r_alpha(1)),
                ((39 - sqrt21) / 120, _r_alpha((1 - sqrt21) / 10)),
            ]
        elif index == 7:
            self.degree = 3
            data = [
                (frac(9, 20), _s3(symbolic)),
                (frac(1, 20), _r_alpha(1)),
                (frac(2, 15), _r_alpha(-frac(1, 2))),
            ]
        elif index == 8:
            self.degree = 4
            sqrt10 = sqrt(10)
            sqrt_b = sqrt(950 - 220 * sqrt10)
            a1 = (-10 + 5 * sqrt10 + sqrt_b) / 30
            a2 = (-10 + 5 * sqrt10 - sqrt_b) / 30
            data = [
                ((5 * a2 - 2) / (60 * a1 ** 2 * (a2 - a1)), _r_alpha(a1)),
                ((5 * a1 - 2) / (60 * a2 ** 2 * (a1 - a2)), _r_alpha(a2)),
            ]
        elif index == 9:
            self.degree = 4
            data = [
                (frac(27, 80), _s3(symbolic)),
                (frac(8, 105), _r_alpha(-frac(1, 2))),
                (frac(81, 560), _r_alpha(frac(2, 3))),
            ]
        elif index == 10:
            self.degree = 4
            sqrt13 = sqrt(13)
            data = [
                ((11 - 1 * sqrt13) / 360, _r_alpha(1)),
                ((80 - 16 * sqrt13) / 360, _r_alpha(-frac(1, 2))),
                ((29 + 17 * sqrt13) / 360, _r_alpha((-1 + sqrt13) / 6)),
            ]
        elif index == 11:
            self.degree = 4
            sqrt3 = sqrt(3)
            data = [
                (+frac(9, 20), _s3(symbolic)),
                (-frac(1, 60), _r_alpha(1)),
                (+frac(1, 10), _r_gamma_delta((3 + sqrt3) / 6, (3 - sqrt3) / 6)),
            ]
        elif index == 12:
            self.degree = 5
            sqrt15 = sqrt(15)
            data = [
                (frac(9, 40), _s3(symbolic)),
                ((155 - sqrt15) / 1200, _r_alpha((1 + sqrt15) / 7)),
                ((155 + sqrt15) / 1200, _r_alpha((1 - sqrt15) / 7)),
            ]
        else:
            assert index == 13
            self.degree = 5
            data = [
                (frac(81, 320), _s3(symbolic)),
                (frac(1, 90), _r_alpha(1)),
                (frac(16, 225), _r_alpha(-frac(1, 2))),
                (frac(2401, 14400), _r_alpha(frac(4, 7))),
            ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return


def _r_alpha(alpha):
    """From the article:

    mu_i = (1 + (n-1) alpha) / n,
    mu_j = (1 - alpha) / n    for j!=i,

    where n is the number of vertices
    """
    a = (1 + 2 * alpha) / 3
    b = (1 - alpha) / 3
    return numpy.array([[a, b, b], [b, a, b], [b, b, a]])


def _r_gamma_delta(gamma, delta):
    """From the article:

    mu_i = (1 + (n-1) gamma - delta) / n,
    mu_j = (1 + (n-1) delta - gamma) / n,
    mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

    where n is the number of vertices
    """
    a = (1 + 2 * gamma - delta) / 3
    b = (1 + 2 * delta - gamma) / 3
    c = (1 - gamma - delta) / 3
    return numpy.array(
        [[a, b, c], [c, a, b], [b, c, a], [a, c, b], [b, a, c], [c, b, a]]
    )
