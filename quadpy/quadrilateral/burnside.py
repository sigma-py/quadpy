# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s


class Burnside(object):
    '''
    W. Burnside,
    An approximate quadrature formula,
    Messenger of Math., v. 37, 1908, pp. 166-167.
    '''
    def __init__(self):
        self.name = 'Burnside'
        self.degree = 5
        self.weights = numpy.concatenate([
            numpy.full(4, 10.0/49.0),
            numpy.full(4, 9.0/196.0),
            ])
        r = numpy.sqrt(7.0 / 15.0)
        s = numpy.sqrt(7.0 / 9.0)
        self.points = numpy.concatenate([
            _symm_r_0(r),
            _symm_s(s)
            ])

        self.weights *= 4.0
        return
