# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _z


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self):
        self.name = 'Tyler'
        self.degree = 5
        self.weights = numpy.concatenate([
            numpy.full(1, -28.0/45.0),
            numpy.full(4, 1.0/36.0),
            numpy.full(4, 1.0/45.0),
            numpy.full(4, 16.0/45.0),
            ])
        self.points = numpy.concatenate([
            _z(),
            _symm_s(1.0),
            _symm_r_0(1.0),
            _symm_r_0(0.5),
            ])

        self.weights *= 4.0
        return
