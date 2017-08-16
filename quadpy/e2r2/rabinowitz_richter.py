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
                (.1237222328857347e+00, _s40(1.538189001320852)),
                (.6544984694978697e-01, _s4(1.224744871391589)),
                (.5935280476180875e+00, _s4(0.4817165220011443)),
                (.1349017971918148e-02, _s8(2.607349811958554, 0.9663217712794149)),
                ]
        else:
            assert False

        self.points, self.weights = untangle(data)
        return
