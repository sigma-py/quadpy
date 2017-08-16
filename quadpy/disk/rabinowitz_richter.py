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
                (.1851958765246450, _s40(.8377170225998396)),
                (.2930225148631698, _s40(.3924393142315810)),
                (.2296152967863584, _s4(.5505609906724360)),
                (.0387822376116376, _s8(.4249164962326038, .9112013890413142)),
                ]
        else:
            assert False

        self.points, self.weights = untangle(data)
        return
