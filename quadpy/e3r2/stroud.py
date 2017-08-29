# -*- coding: utf-8 -*-
#
from .stroud_secrest import StroudSecrest


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == 'E3r2 5-1':
            self.set_data(StroudSecrest('VII'))
        elif index == 'E3r2 5-2a':
            self.set_data(StroudSecrest('VIIIa'))
        elif index == 'E3r2 5-2b':
            self.set_data(StroudSecrest('VIIIb'))
        elif index == 'E3r2 5-3':
            self.set_data(StroudSecrest('IX'))
        elif index == 'E3r2 7-1a':
            self.set_data(StroudSecrest('Xa'))
        elif index == 'E3r2 7-1b':
            self.set_data(StroudSecrest('Xb'))
        elif index == 'E3r2 7-2a':
            self.set_data(StroudSecrest('XIa'))
        elif index == 'E3r2 7-2b':
            self.set_data(StroudSecrest('XIb'))
        else:
            assert False
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
