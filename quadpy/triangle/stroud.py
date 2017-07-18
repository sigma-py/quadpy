# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s3, _s21, _s111


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, index):
        self.name = 'Stroud({})'.format(index)
        if index == 0:
            self.degree = 1
            self.weights = numpy.array([0.5])
            self.bary = _s3()
        elif index == 1:
            self.degree = 1
            self.weights = numpy.full(3, 1.0/6.0)
            self.bary = _s21(0.0)
        elif index == 2:
            self.degree = 2
            self.weights = numpy.full(3, 1.0/6.0)
            self.bary = _s21(0.5)
        elif index == 3:
            self.degree = 2
            self.weights = numpy.full(3, 1.0/6.0)
            self.bary = _s21(1.0/6.0)
        elif index == 4:
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(1, -0.28125),
                numpy.full(3, 25.0/96.0)
                ])
            self.bary = numpy.concatenate([
                _s3(),
                _s21(0.2)
                ])
        elif index == 5:
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(3, -0.0149273987289826666941168633072000),
                numpy.full(3, 0.181594065395649333360783529973866)
                ])
            self.bary = numpy.concatenate([
                _s21(0.0),
                _s21(0.147247476834805333113731760209066)
                ])
        elif index == 6:
            self.degree = 3
            self.weights = numpy.full(6, 1.0/12.0)
            self.bary = _s111(
                0.109039009072877212324834667569910,
                0.231933368553030572496784561174692
                )
        elif index == 7:
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(3, 1.0/60.0),
                numpy.full(3, 0.15)
                ])
            self.bary = numpy.concatenate([
                _s21(0.5),
                _s21(1.0/6.0)
                ])
        elif index == 8:
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(1, 0.225),
                numpy.full(3, 1.0/15.0),
                numpy.full(3, 0.025)
                ])
            self.bary = numpy.concatenate([
                _s3(),
                _s21(0.5),
                _s21(0.0)
                ])
        else:
            assert index == 9
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(1, 0.1125),
                numpy.full(3, 0.0629695902724135762978419727500906),
                numpy.full(3, 0.0661970763942530903688246939165759)
                ])
            self.bary = numpy.concatenate([
                _s3(),
                _s21(0.101286507323456338800987361915123),
                _s21(0.470142064105115089770441209513447)
                ])

        self.weights *= 2
        self.points = self.bary[:, 1:]
        return
