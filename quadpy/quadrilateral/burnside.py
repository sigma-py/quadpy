# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s

from ..helpers import untangle


class Burnside(object):
    '''
    W. Burnside,
    An approximate quadrature formula,
    Messenger of Math., v. 37, 1908, pp. 166-167.
    '''
    def __init__(self):
        self.name = 'Burnside'
        self.degree = 5
        r = numpy.sqrt(7.0 / 15.0)
        s = numpy.sqrt(7.0 / 9.0)
        data = [
            (10.0/49.0, _symm_r_0(r)),
            (9.0/196.0, _symm_s(s)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0
        return
