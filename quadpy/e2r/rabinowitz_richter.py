# -*- coding: utf-8 -*-
#
from .helpers import _s40, _s8, _s4
from ..helpers import untangle


class RabinowitzRichter(object):
    '''
    Philip Rabinowitz and Nira Richter,
    Perfectly Symmetric Two-Dimensional Integration Formulas with Minimal
    Numbers of Points,
    Mathematics of Computation, Vol. 23, No. 108 (Oct., 1969), pp. 765-779,
    <https://dx.doi.org/10.2307/2004962>.
    '''
    def __init__(self, index):
        self.name = 'RabinowitzRichter({})'.format(index)
        if index == 1:
            self.degree = 9
            data = [
                (.3380228176732269e-1, _s40(6.822859174233539)),
                (.1467201651910359e+1, _s40(1.901350903458987)),
                (.6973178170307865e-1, _s4(4.260195453867070)),
                (.3030570706813315e-4, _s8(6.693991707281686, 14.77112509749386)),
                ]
        else:
            assert False

        self.points, self.weights = untangle(data)
        return
