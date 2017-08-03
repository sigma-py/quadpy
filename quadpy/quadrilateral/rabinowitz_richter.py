# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _symm_s_t


class RabinowitzRichter(object):
    '''
    Philip Rabinowitz and Nira Richter,
    Perfectly Symmetric Two-Dimensional Integration Formulas with Minimal
    Numbers of Points,
    Mathematics of Computation, Vol. 23, No. 108 (Oct., 1969), pp. 765-779,
    <https://dx.doi.org/10.2307/2004962>.
    '''
    def __init__(self, index):
        self.name = 'RabinowitzRichter({})'.format(index)
        if index == 1:
            self.degree = 9
            r1 = 0.9845398119422523
            r2 = 0.4888863428423724
            r3 = 0.9395672874215217
            r4 = 0.8367103250239890

            s4 = 0.5073767736746132

            B1 = 0.0716134247098111
            B2 = 0.4540903525515453
            B3 = 0.0427846154667780
            B4 = 0.2157558036359328

            # self.weights = numpy.concatenate([
            #     numpy.full(4, B1),
            #     numpy.full(4, B2),
            #     numpy.full(4, B3),
            #     numpy.full(8, B4),
            #     ])
            self.weights = numpy.repeat([B1, B2, B3, B4], [4, 4, 4, 8])
            self.points = numpy.concatenate([
                _symm_r_0(r1),
                _symm_r_0(r2),
                _symm_s(r3),
                _symm_s_t(r4, s4),
                ])
        else:
            assert False

        return
