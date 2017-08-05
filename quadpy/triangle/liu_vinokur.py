# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s3

from ..helpers import untangle


class LiuVinokur(object):
    '''
    Y. Liu and M. Vinokur,
    Exact Integrations of Polynomials and Symmetric Quadrature Formulas over
    Arbitrary Polyhedral Grids,
    Journal of Computational Physics, 140, 122–147 (1998).
    DOI: 10.1006/jcph.1998.5884,
    <https://dx.doi.org/10.1006/jcph.1998.5884>.
    '''
    def __init__(self, index):
        self.name = 'LV(%d)' % index
        if index == 1:
            self.degree = 1
            data = [(1.0, _s3())]
        elif index == 2:
            self.degree = 1
            data = [(1.0/3.0, _r_alpha(1.0))]
        elif index == 3:
            self.degree = 2
            data = [
                (1.0/3.0, _r_alpha(-0.5)),
                ]
        elif index == 4:
            self.degree = 2
            data = [
                (0.75, _s3()),
                (1.0/12.0, _r_alpha(1.0)),
                ]
        elif index == 5:
            self.degree = 3
            data = [
                (-9.0/16.0, _s3()),
                # Wrongly specified in the article as 25 (instead of 2/5).
                (25.0/48.0, _r_alpha(0.4)),
                ]
        elif index == 6:
            self.degree = 3
            sqrt21 = numpy.sqrt(21.0)
            data = [
                ((1.00 + sqrt21)/120.0, _r_alpha(1.0)),
                ((39.0 - sqrt21)/120.0, _r_alpha((1.0 - sqrt21) / 10.0)),
                ]
        elif index == 7:
            self.degree = 3
            data = [
                (9.0/20.0, _s3()),
                (1.0/20.0, _r_alpha(1.0)),
                (2.0/15.0, _r_alpha(-0.5)),
                ]
        elif index == 8:
            self.degree = 4
            sqrt10 = numpy.sqrt(10.0)
            sqrt_b = numpy.sqrt(950.0 - 220.0*sqrt10)
            a1 = (-10.0 + 5.0*sqrt10 + sqrt_b) / 30.0
            a2 = (-10.0 + 5.0*sqrt10 - sqrt_b) / 30.0
            data = [
                ((5*a2-2) / (60*a1**2 * (a2 - a1)), _r_alpha(a1)),
                ((5*a1-2) / (60*a2**2 * (a1 - a2)), _r_alpha(a2)),
                ]
        elif index == 9:
            self.degree = 4
            data = [
                (27.0/80.00, _s3()),
                (8.00/105.0, _r_alpha(-0.5)),
                (81.0/560.0, _r_alpha(2.0/3.0)),
                ]
        elif index == 10:
            self.degree = 4
            sqrt13 = numpy.sqrt(13.0)
            data = [
                ((11.0 - 1.00*sqrt13)/360.0, _r_alpha(1.0)),
                ((80.0 - 16.0*sqrt13)/360.0, _r_alpha(-0.5)),
                ((29.0 + 17.0*sqrt13)/360.0, _r_alpha((-1.0 + sqrt13) / 6.0)),
                ]
        elif index == 11:
            self.degree = 4
            sqrt3 = numpy.sqrt(3.0)
            data = [
                (0.45, _s3()),
                (-1.0/60.0, _r_alpha(1.0)),
                (0.1, _r_gamma_delta((3.0 + sqrt3)/6.0, (3.0 - sqrt3)/6.0)),
                ]
        elif index == 12:
            self.degree = 5
            sqrt15 = numpy.sqrt(15.0)
            data = [
                (9.0/40.0, _s3()),
                ((155.0 - sqrt15)/1200.0, _r_alpha((1.0 + sqrt15) / 7.0)),
                ((155.0 + sqrt15)/1200.0, _r_alpha((1.0 - sqrt15) / 7.0)),
                ]
        else:
            assert index == 13
            self.degree = 5
            data = [
                (81.0/320.0, _s3()),
                (1.0/90.0, _r_alpha(1.0)),
                (16.0/225.0, _r_alpha(-0.5)),
                (2401.0/14400.0, _r_alpha(4.0/7.0)),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return


def _r_alpha(alpha):
    '''From the article:

    mu_i = (1 + (n-1) alpha) / n,
    mu_j = (1 - alpha) / n    for j!=i,

    where n is the number of vertices
    '''
    a = (1.0 + 2*alpha) / 3.0
    b = (1.0 - alpha) / 3.0
    return numpy.array([
        [a, b, b],
        [b, a, b],
        [b, b, a],
        ])


def _r_gamma_delta(gamma, delta):
    '''From the article:

    mu_i = (1 + (n-1) gamma - delta) / n,
    mu_j = (1 + (n-1) delta - gamma) / n,
    mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

    where n is the number of vertices
    '''
    a = (1.0 + 2*gamma - delta) / 3.0
    b = (1.0 + 2*delta - gamma) / 3.0
    c = (1.0 - gamma - delta) / 3.0
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [a, c, b],
        [b, a, c],
        [c, b, a],
        ])
