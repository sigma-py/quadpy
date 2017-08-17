# -*- coding: utf-8 -*-
#
from ..helpers import untangle, fsd
from .helpers import integrate_monomial_over_unit_nsphere


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, n, index):
        self.dim = n
        if index == 'Un 3-1':
            self.degree = 3
            data = [
                (0.5/n, fsd(n, 1.0, 1)),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
