# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
from ..helpers import untangle


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
        self.degree = 8
        data = [
            (0.0378610912003147, _s3()),
            (0.0376204254131829, _s21(0.1012865073234563)),
            (0.0783573522441174, _s21(0.4701420641051151)),
            (0.1162714796569659, _s21(0.2321023267750504)),
            (0.0134442673751655, _s21(0.0294808608844396)),
            (0.0375097224552317, _s111(0.7384168123405100, 0.2321023267750504)),
            ]
        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
