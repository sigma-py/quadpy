# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s31, _s22


class Gatermann(object):
    '''
    Karin Gatermann,
    Linear Representations of Finite Groups and The Ideal Theoretical
    Construction of G-Invariant Cubature Formulas,
    Numerical Integration pp 25-35,
    Part of the NATO ASI Series book series (ASIC, volume 357).
    '''
    def __init__(self):
        self.name = 'Gatermann'
        self.degree = 5
        self.weights = 6 * numpy.concatenate([
            numpy.full(4, 9.73033316198362119165356216965707e-06),
            numpy.full(4, 8.99031481668747219698547129902142e-03),
            numpy.full(6, 2.17777476778781405656596945369837e-02),
            ])
        bary = numpy.concatenate([
            _s31(0.656936552995394536166881327385593),
            _s31(0.0801424420792727848879183805550907),
            _s22(0.404475329343454044779549906725159),
            ])

        self.points = bary[:, 1:]
        return
