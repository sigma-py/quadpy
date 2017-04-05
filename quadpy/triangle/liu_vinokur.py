# -*- coding: utf-8 -*-
#
from .helpers import _s3
import numpy


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
            self.weights = numpy.concatenate([
                1.0 * numpy.ones(1),
                ])
            bary = numpy.concatenate([
                _s3(),
                ])
            self.degree = 1
        elif index == 2:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(1.0),
                ])
            self.degree = 1
        elif index == 3:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(-0.5),
                ])
            self.degree = 2
        elif index == 4:
            self.weights = numpy.concatenate([
                0.75 * numpy.ones(1),
                1.0/12.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                ])
            self.degree = 2
        elif index == 5:
            self.weights = numpy.concatenate([
                -9.0/16.0 * numpy.ones(1),
                25.0/48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                # Wrongly specified in the article as 25 (instead of 2/5).
                self._r_alpha(0.4),
                ])
            self.degree = 3
        elif index == 6:
            self.weights = numpy.concatenate([
                (1.0 + numpy.sqrt(21.0)) / 120.0 * numpy.ones(3),
                (39.0 - numpy.sqrt(21.0)) / 120.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(1.0),
                self._r_alpha((1.0 - numpy.sqrt(21.0)) / 10.0),
                ])
            self.degree = 3
        elif index == 7:
            self.weights = numpy.concatenate([
                9.0/20.0 * numpy.ones(1),
                1.0/20.0 * numpy.ones(3),
                2.0/15.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                self._r_alpha(-0.5),
                ])
            self.degree = 3
        elif index == 8:
            sqrt10 = numpy.sqrt(10)
            alpha1 = (-10 + 5*sqrt10 + numpy.sqrt(950.0 - 220*sqrt10)) / 30.0
            alpha2 = (-10 + 5*sqrt10 - numpy.sqrt(950.0 - 220*sqrt10)) / 30.0
            self.weights = numpy.concatenate([
                (5*alpha2-2) / (60*alpha1**2 * (alpha2 - alpha1))
                * numpy.ones(3),
                (5*alpha1-2) / (60*alpha2**2 * (alpha1 - alpha2))
                * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(alpha1),
                self._r_alpha(alpha2),
                ])
            self.degree = 4
        elif index == 9:
            self.weights = numpy.concatenate([
                27.0/80.0 * numpy.ones(1),
                8.0/105.0 * numpy.ones(3),
                81.0/560.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(-0.5),
                self._r_alpha(2.0/3.0),
                ])
            self.degree = 4
        elif index == 10:
            self.weights = numpy.concatenate([
                (11.0 - numpy.sqrt(13)) / 360.0 * numpy.ones(3),
                (80.0 - 16*numpy.sqrt(13)) / 360.0 * numpy.ones(3),
                (29.0 + 17*numpy.sqrt(13)) / 360.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(1.0),
                self._r_alpha(-0.5),
                self._r_alpha((-1.0 + numpy.sqrt(13.0)) / 6.0),
                ])
            self.degree = 4
        elif index == 11:
            self.weights = numpy.concatenate([
                0.45 * numpy.ones(1),
                -1.0 / 60.0 * numpy.ones(3),
                0.1 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                self._r_gamma_delta(
                    (3.0 + numpy.sqrt(3.0)) / 6.0,
                    (3.0 - numpy.sqrt(3.0)) / 6.0
                    ),
                ])
            self.degree = 4
        elif index == 12:
            self.weights = numpy.concatenate([
                9.0/40.0 * numpy.ones(1),
                (155.0 - numpy.sqrt(15.0))/1200.0 * numpy.ones(3),
                (155.0 + numpy.sqrt(15.0))/1200.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha((1.0 + numpy.sqrt(15.0)) / 7.0),
                self._r_alpha((1.0 - numpy.sqrt(15.0)) / 7.0),
                ])
            self.degree = 5
        else:
            assert index == 13
            self.weights = numpy.concatenate([
                81.0/320.0 * numpy.ones(1),
                1.0/90.0 * numpy.ones(3),
                16.0/225.0 * numpy.ones(3),
                2401.0/14400.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                self._r_alpha(-0.5),
                self._r_alpha(4.0/7.0),
                ])
            self.degree = 5

        self.points = bary[:, 1:]
        return

    def _r_alpha(self, alpha):
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

    def _r_gamma_delta(self, gamma, delta):
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
