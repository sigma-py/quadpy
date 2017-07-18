# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s3, _s21, _s111


class Cubtri(object):
    '''
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    Laurie, D. P.,
    Algorithm 584: CUBTRI: Automatic Cubature over a Triangle,
    ACM Trans. Math. Softw.,
    June 1982,
    <http://dl.acm.org/citation.cfm?id=356001>.
    '''
    def __init__(self):
        self.name = 'CUBTRI'
        self.weights = numpy.concatenate([
            numpy.full(1, 0.0378610912003147),
            numpy.full(3, 0.0376204254131829),
            numpy.full(3, 0.0783573522441174),
            numpy.full(3, 0.1162714796569659),
            numpy.full(3, 0.0134442673751655),
            numpy.full(6, 0.0375097224552317),
            ])

        self.bary = numpy.concatenate([
            _s3(),
            _s21(0.1012865073234563),
            _s21(0.4701420641051151),
            _s21(0.2321023267750504),
            _s21(0.0294808608844396),
            _s111(0.7384168123405100, 0.2321023267750504),
            ])
        self.points = self.bary[:, 1:]
        self.degree = 8
        return
