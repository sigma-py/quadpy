# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s4, _s31


class Zienkiewicz(object):
    '''
    Olgierd Zienkiewicz,
    The Finite Element Method,
    Sixth Edition,
    Butterworth-Heinemann, 2005,
    ISBN: 0750663200,
    <http://www.sciencedirect.com/science/book/9780750664318>,
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html>.
    '''
    def __init__(self, index):
        if index == 4:
            self.weights = 0.25 * numpy.ones(4)
            bary = _s31(0.1381966011250105)
            self.degree = 2
        else:
            assert index == 5
            self.weights = numpy.concatenate([
                -0.8 * numpy.ones(1),
                0.45 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/6.0),
                ])
            self.degree = 3

        self.points = bary[:, [1, 2, 3]]
        return
