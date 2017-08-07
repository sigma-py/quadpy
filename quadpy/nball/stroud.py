# -*- coding: utf-8 -*-
#
import numpy

from . import hammer_stroud
from . import stroud1957

from .helpers import volume_unit_ball
from ..helpers import pm, untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index):
        self.name = 'Stroud({})'.format(index)
        self.dim = n
        if index == 'Sn 2-1':
            self.set_data(stroud1957.Stroud1957(n))
        elif index == 'Sn 3-1':
            self.set_data(hammer_stroud.HammerStroud(n, alpha=0.0))
        elif index == 'Sn 3-2':
            self.degree = 3
            r = numpy.sqrt(1.0 / (n+2))
            data = [
                (0.5**n, pm(n, r))
                ]
            self.points, self.weights = untangle(data)
            self.weights *= volume_unit_ball(n)
        else:
            assert False

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
