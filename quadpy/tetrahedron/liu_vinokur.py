# -*- coding: utf-8 -*-
#
import numpy
from sympy import sqrt, Rational as fr

from .helpers import _s4

# TODO update weight/points specification
class LiuVinokur(object):
    '''
    Y. Liu and M. Vinokur,
    Exact Integrations of Polynomials and Symmetric Quadrature Formulas over
    Arbitrary Polyhedral Grids,
    Journal of Computational Physics, 140, 122–147 (1998).
    DOI: 10.1006/jcph.1998.5884,
    <http://dx.doi.org/10.1006/jcph.1998.5884>.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.concatenate([
                numpy.full(1, 1),
                ])
            bary = numpy.concatenate([
                _s4(),
                ])
            self.degree = 1
        elif index == 2:
            self.weights = numpy.concatenate([
                numpy.full(4, fr(1, 4)),
                ])
            bary = numpy.concatenate([
                _r_alpha(1.0),
                ])
            self.degree = 1
        elif index == 3:
            self.weights = numpy.concatenate([
                numpy.full(4, fr(1, 4)),
                ])
            bary = numpy.concatenate([
                _r_alpha(1/sqrt(5)),
                ])
            self.degree = 2
        elif index == 4:
            self.weights = numpy.concatenate([
                numpy.full(1, fr(4, 5)),
                numpy.full(4, fr(1, 20)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1),
                ])
            self.degree = 2
        elif index == 5:
            self.weights = numpy.concatenate([
                numpy.full(1, -fr(4, 5)),
                numpy.full(4, fr(9, 20)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(fr(1, 3)),
                ])
            self.degree = 3
        elif index == 6:
            self.weights = numpy.concatenate([
                numpy.full(4, fr(1, 40)),
                numpy.full(4, fr(9, 40)),
                ])
            bary = numpy.concatenate([
                _r_alpha(1),
                _r_alpha(-fr(1, 3)),
                ])
            self.degree = 3
        elif index == 7:
            self.weights = numpy.concatenate([
                numpy.full(1, -fr(148, 1875)),
                numpy.full(4, fr(343, 7500)),
                numpy.full(6, fr(56, 375)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(fr(5, 7)),
                _r_beta(sqrt(70)/28),
                ])
            self.degree = 4
        elif index == 8:
            alpha1 = (
                + sqrt(65944 - 19446*sqrt(11))
                + 51*sqrt(11) - 154
                ) / 89
            alpha2 = (
                - sqrt(65944 - 19446*sqrt(11))
                + 51*sqrt(11) - 154
                ) / 89
            self.weights = numpy.concatenate([
                numpy.full(
                    4, (17*alpha2 - 7)/(420*alpha1**2 * (alpha2 - alpha1))
                    ),
                numpy.full(
                    4, (17*alpha1 - 7)/(420*alpha2**2 * (alpha1 - alpha2))
                    ),
                numpy.full(6, fr(2, 105)),
                ])
            bary = numpy.concatenate([
                _r_alpha(alpha1),
                _r_alpha(alpha2),
                _r_beta(fr(1, 2)),
                ])
            self.degree = 4
        elif index == 9:
            self.weights = numpy.concatenate([
                numpy.full(1, -fr(32, 15)),
                numpy.full(4, fr(3, 280)),
                numpy.full(4, fr(125, 168)),
                numpy.full(6, fr(2, 105)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1),
                _r_alpha(fr(1, 5)),
                _r_beta(fr(1, 2)),
                ])
            self.degree = 4
        elif index == 10:
            self.weights = numpy.concatenate([
                numpy.full(1, fr(32, 105)),
                numpy.full(4, -fr(31, 840)),
                numpy.full(4, fr(27, 280)),
                numpy.full(12, fr(4, 105)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1),
                _r_alpha(-fr(1, 3)),
                _r_gamma_delta(
                    (2 + sqrt(2)) / 4,
                    (2 - sqrt(2)) / 4,
                    ),
                ])
            self.degree = 4
        elif index == 11:
            self.weights = numpy.concatenate([
                (11 - 4*sqrt(2)) / numpy.full(4, 840),
                (243 - 108*sqrt(2)) / numpy.full(4, 1960),
                (62 + 44*sqrt(2)) / numpy.full(4, 735),
                numpy.full(6, fr(2, 105)),
                ])
            bary = numpy.concatenate([
                _r_alpha(1),
                _r_alpha(-fr(1, 3)),
                _r_alpha(sqrt(2) - 1),
                _r_beta(fr(1, 2)),
                ])
            self.degree = 4
        elif index == 12:
            lmbda = fr(4, 27) * (4 * sqrt(79)*numpy.cos(
                (numpy.arccos(67*sqrt(79)/24964) + 2*numpy.pi) / 3
                ) + 71
                )
            alpha1 = (
                + sqrt(9*lmbda**2 - 248*lmbda + 1680) + 28 - 3*lmbda
                ) / (112 - 10*lmbda)
            alpha2 = (
                - sqrt(9*lmbda**2 - 248*lmbda + 1680) + 28 - 3*lmbda
                ) / (112 - 10*lmbda)
            w1 = ((21 - lmbda)*alpha2 - 7) \
                / (420*alpha1**2 * (alpha2 - alpha1))
            w2 = ((21 - lmbda)*alpha1 - 7) \
                / (420*alpha2**2 * (alpha1 - alpha2))
            self.weights = numpy.concatenate([
                numpy.full(4, w1),
                numpy.full(4, w2),
                numpy.full(6, lmbda**2 / 840),
                ])
            bary = numpy.concatenate([
                _r_alpha(alpha1),
                _r_alpha(alpha2),
                _r_beta(1/sqrt(lmbda)),
                ])
            self.degree = 5
        elif index == 13:
            self.weights = numpy.concatenate([
                numpy.full(1, -fr(16, 21)),
                numpy.full(4, (2249 - 391*sqrt(13)) / 10920),
                numpy.full(4, (2249 + 391*sqrt(13)) / 10920),
                numpy.full(6, fr(2, 105)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha((2 + sqrt(13)) / 9),
                _r_alpha((2 - sqrt(13)) / 9),
                _r_beta(fr(1, 2)),
                ])
            self.degree = 5
        else:
            assert index == 14
            self.weights = numpy.concatenate([
                numpy.full(1, fr(16, 105)),
                numpy.full(4, fr(1, 280)),
                numpy.full(4, fr(81, 1400)),
                numpy.full(4, fr(64, 525)),
                numpy.full(6, fr(2, 105)),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1),
                _r_alpha(-fr(1, 3)),
                _r_alpha(fr(1, 2)),
                _r_beta(fr(1, 2)),
                ])
            self.degree = 5

        self.points = bary[:, 1:]
        return


def _r_alpha(alpha):
    '''From the article:

    mu_i = (1 + (n-1) alpha) / n,
    mu_j = (1 - alpha) / n    for j!=i,

    where n is the number of vertices
    '''
    a = (1 + 3*alpha) / 4
    b = (1 - alpha) / 4
    return numpy.array([
        [a, b, b, b],
        [b, a, b, b],
        [b, b, a, b],
        [b, b, b, a],
        ])


def _r_beta(beta):
    '''From the article:

    mu_i = (1+(n-2)*beta) / n,
    mu_j = mu_i,
    mu_k = (1 - 2*beta) / n    for k!=i, k!=j,

    where n is the number of vertices.
    '''
    a = (1 + 2*beta) / 4
    b = (1 - 2*beta) / 4
    return numpy.array([
        [a, a, b, b],
        [a, b, a, b],
        [b, a, a, b],
        [a, b, b, a],
        [b, a, b, a],
        [b, b, a, a],
        ])


def _r_gamma_delta(gamma, delta):
    '''From the article:

    mu_i = (1 + (n-1) gamma - delta) / n,
    mu_j = (1 + (n-1) delta - gamma) / n,
    mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

    where n is the number of vertices
    '''
    b = (1 + 3*gamma - delta) / 4
    c = (1 + 3*delta - gamma) / 4
    a = (1 - gamma - delta) / 4
    return numpy.array([
        [a, a, b, c],
        [a, b, a, c],
        [b, a, a, c],
        [a, b, c, a],
        [b, a, c, a],
        [b, c, a, a],
        [a, a, c, b],
        [a, c, a, b],
        [c, a, a, b],
        [a, c, b, a],
        [c, a, b, a],
        [c, b, a, a],
        ])
