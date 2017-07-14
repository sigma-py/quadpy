# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s3, _s21


class SevenPoint(object):
    def __init__(self):
        self.weights = numpy.concatenate([
            numpy.full(1, 0.45),
            numpy.full(3, 0.05),
            2.0 / numpy.full(3, 15.0),
            ])
        bary = numpy.concatenate([
            _s3(),
            _s21(0.0),
            _s21(0.5),
            ])
        self.points = bary[:, [1, 2]]
        self.degree = 3
        self.name = 'seven-point'
        return
