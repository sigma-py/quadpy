# -*- coding: utf-8 -*-
#
import numpy

from .mustard_lyness_blatt import MustardLynessBlatt
from .tyler import Tyler

from .helpers import fs_rrr

from ..ncube import Ewing
from ..helpers import untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, index):
        reference_volume = 8.0
        if index == 'C3 3-1':
            self.set_data(Tyler())
        elif index == 'C3 3-2':
            # product Gauss
            self.degree = 3
            data = [
                (1.0/8.0, fs_rrr(numpy.sqrt(1.0/3.0)))
                ]
            self.points, self.weights = untangle(data)
            self.weights *= reference_volume
        elif index == 'C3 3-3':
            self.set_data(Ewing(3))
        elif index == 'C3 3-4':
            self.set_data(MustardLynessBlatt())
        else:
            assert False

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
