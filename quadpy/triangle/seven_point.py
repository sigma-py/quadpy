# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21
import numpy


class SevenPoint(object):
    def __init__(self):
        self.weights = numpy.concatenate([
            0.45 * numpy.ones(1),
            0.05 * numpy.ones(3),
            2.0 / 15.0 * numpy.ones(3),
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
