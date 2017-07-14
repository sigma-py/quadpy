# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s4


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
                numpy.full(1, 1.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                ])
            self.degree = 1
        elif index == 2:
            self.weights = numpy.concatenate([
                numpy.full(4, 0.25),
                ])
            bary = numpy.concatenate([
                _r_alpha(1.0),
                ])
            self.degree = 1
        elif index == 3:
            self.weights = numpy.concatenate([
                numpy.full(4, 0.25),
                ])
            bary = numpy.concatenate([
                _r_alpha(1.0 / numpy.sqrt(5.0)),
                ])
            self.degree = 2
        elif index == 4:
            self.weights = numpy.concatenate([
                numpy.full(1, 0.8),
                numpy.full(4, 0.05),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1.0),
                ])
            self.degree = 2
        elif index == 5:
            self.weights = numpy.concatenate([
                numpy.full(1, -0.8),
                numpy.full(4, 0.45),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1.0/3.0),
                ])
            self.degree = 3
        elif index == 6:
            self.weights = numpy.concatenate([
                numpy.full(4, 1.0/40.0),
                numpy.full(4, 9.0/40.0),
                ])
            bary = numpy.concatenate([
                _r_alpha(1.0),
                _r_alpha(-1.0/3.0),
                ])
            self.degree = 3
        elif index == 7:
            self.weights = numpy.concatenate([
                numpy.full(1, -148.0/1875.0),
                numpy.full(4, 343.0/7500.0),
                numpy.full(6, 56.0/375.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(5.0/7.0),
                _r_beta(numpy.sqrt(70.0)/28.0),
                ])
            self.degree = 4
        elif index == 8:
            alpha1 = (
                + numpy.sqrt(65944.0 - 19446*numpy.sqrt(11))
                + 51*numpy.sqrt(11) - 154.0
                ) / 89.0
            alpha2 = (
                - numpy.sqrt(65944.0 - 19446*numpy.sqrt(11))
                + 51*numpy.sqrt(11) - 154.0
                ) / 89.0
            self.weights = numpy.concatenate([
                numpy.full(
                    4, (17*alpha2 - 7.0)/(420.0*alpha1**2 * (alpha2 - alpha1))
                    ),
                numpy.full(
                    4, (17*alpha1 - 7.0)/(420.0*alpha2**2 * (alpha1 - alpha2))
                    ),
                numpy.full(6, 2.0/105.0),
                ])
            bary = numpy.concatenate([
                _r_alpha(alpha1),
                _r_alpha(alpha2),
                _r_beta(0.5),
                ])
            self.degree = 4
        elif index == 9:
            self.weights = numpy.concatenate([
                numpy.full(1, -32.0/15.0),
                numpy.full(4, 3.0/280.0),
                numpy.full(4, 125.0/168.0),
                numpy.full(6, 2.0/105.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1),
                _r_alpha(0.2),
                _r_beta(0.5),
                ])
            self.degree = 4
        elif index == 10:
            self.weights = numpy.concatenate([
                numpy.full(1, 32.0/105.0),
                numpy.full(4, -31.0/840.0),
                numpy.full(4, 27.0/280.0),
                numpy.full(12, 4.0/105.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1),
                _r_alpha(-1.0/3.0),
                _r_gamma_delta(
                    (2 + numpy.sqrt(2.0)) / 4.0,
                    (2 - numpy.sqrt(2.0)) / 4.0,
                    ),
                ])
            self.degree = 4
        elif index == 11:
            self.weights = numpy.concatenate([
                (11.0 - 4*numpy.sqrt(2.0)) / numpy.full(4, 840.0),
                (243.0 - 108*numpy.sqrt(2.0)) / numpy.full(4, 1960.0),
                (62.0 + 44*numpy.sqrt(2.0)) / numpy.full(4, 735.0),
                numpy.full(6, 2.0/105.0),
                ])
            bary = numpy.concatenate([
                _r_alpha(1),
                _r_alpha(-1.0/3.0),
                _r_alpha(numpy.sqrt(2.0) - 1.0),
                _r_beta(0.5),
                ])
            self.degree = 4
        elif index == 12:
            lmbda = 4.0/27.0 * (4.0 * numpy.sqrt(79.0)*numpy.cos(
                (numpy.arccos(67*numpy.sqrt(79.0)/24964.0) + 2*numpy.pi) / 3.0
                ) + 71.0
                )
            alpha1 = (
                + numpy.sqrt(9*lmbda**2 - 248*lmbda + 1680) + 28.0 - 3*lmbda
                ) / (112.0 - 10*lmbda)
            alpha2 = (
                - numpy.sqrt(9*lmbda**2 - 248*lmbda + 1680) + 28.0 - 3*lmbda
                ) / (112.0 - 10*lmbda)
            w1 = ((21.0 - lmbda)*alpha2 - 7.0) \
                / (420.0*alpha1**2 * (alpha2 - alpha1))
            w2 = ((21.0 - lmbda)*alpha1 - 7.0) \
                / (420.0*alpha2**2 * (alpha1 - alpha2))
            self.weights = numpy.concatenate([
                numpy.full(4, w1),
                numpy.full(4, w2),
                numpy.full(6, lmbda**2 / 840.0),
                ])
            bary = numpy.concatenate([
                _r_alpha(alpha1),
                _r_alpha(alpha2),
                _r_beta(1.0 / numpy.sqrt(lmbda)),
                ])
            self.degree = 5
        elif index == 13:
            self.weights = numpy.concatenate([
                numpy.full(1, -16.0/21.0),
                (2249.0 - 391.0*numpy.sqrt(13.0)) / numpy.full(4, 10920.0),
                (2249.0 + 391.0*numpy.sqrt(13.0)) / numpy.full(4, 10920.0),
                2.0 / numpy.full(6, 105.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha((2.0 + numpy.sqrt(13.0)) / 9.0),
                _r_alpha((2.0 - numpy.sqrt(13.0)) / 9.0),
                _r_beta(0.5),
                ])
            self.degree = 5
        else:
            assert index == 14
            self.weights = numpy.concatenate([
                numpy.full(1, 16.0/105.0),
                numpy.full(4, 1.0/280.0),
                numpy.full(4, 81.0/1400.0),
                numpy.full(4, 64.0/525.0),
                numpy.full(6, 2.0/105.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _r_alpha(1.0),
                _r_alpha(-1.0/3.0),
                _r_alpha(0.5),
                _r_beta(0.5),
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
    a = (1.0 + 3*alpha) / 4.0
    b = (1.0 - alpha) / 4.0
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
    a = (1.0 + 2*beta) / 4.0
    b = (1.0 - 2*beta) / 4.0
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
    b = (1.0 + 3*gamma - delta) / 4.0
    c = (1.0 + 3*delta - gamma) / 4.0
    a = (1.0 - gamma - delta) / 4.0
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
