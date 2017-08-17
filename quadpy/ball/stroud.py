# -*- coding: utf-8 -*-
#
from .ditkin import Ditkin
from .hammer_stroud import HammerStroud


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == 'S3 3-1':
            self.set_data(HammerStroud('11-3'))
        elif index == 'S3 5-1':
            self.set_data(Ditkin(1))
        elif index == 'S3 5-2':
            self.set_data(Ditkin(2))
        elif index == 'S3 7-3':
            self.set_data(Ditkin(3))
        else:
            assert False

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
