# -*- coding: utf-8 -*-
#
from .tyler import Tyler


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.

    <https://people.sc.fsu.edu/~jburkardt/m_src/stroud/square_unit_set.m>
    <http://nines.cs.kuleuven.be/ecf/mtables.html>
    '''
    def __init__(self, index):
        if index == 'C3 3-1':
            self.set_data(Tyler())
        else:
            assert False

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
