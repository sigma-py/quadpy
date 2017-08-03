# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _z, _symm_s_t


class Meister(object):
    '''
    Bernd Meister,
    On a Family of Cubature Formulae,
    Comput J (1966) 8 (4): 368-371,
    <https://doi.org/10.1093/comjnl/8.4.368>.
    '''
    def __init__(self):
        self.name = 'Meister'
        self.degree = 7

        r = 2.0/3.0
        s = 1.0/3.0

        self.weights = numpy.concatenate([
            numpy.full(1, 1024.0/6720.0),
            numpy.full(4, 576.0/6720.0),
            numpy.full(4, 576.0/6720.0),
            numpy.full(4, -9.0/6720.0),
            numpy.full(8, 117.0/6720.0),
            numpy.full(4, 47.0/6720.0),
            ])
        self.points = numpy.concatenate([
            _z(),
            _symm_s(r),
            _symm_r_0(r),
            _symm_s(s),
            _symm_s_t(1.0, s),
            _symm_s(1.0),
            ])

        self.weights *= 4.0
        return
