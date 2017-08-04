# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _pm2

from ..helpers import untangle


class Phillips(object):
    '''
    G.M. Phillips,
    Numerical integration in two and three dimensions,
    Comput J (1967) 10 (2): 202-204,
    <https://doi.org/10.1093/comjnl/10.2.202>.

    Abtract:
    Gaussian-type quadrature formulae are derived for a rectangular region of
    two or three dimensions.
    '''
    def __init__(self):
        self.name = 'Phillips'

        c = 3.0*numpy.sqrt(385.0)
        r = numpy.sqrt((105.0 + c) / 140.0)
        s = numpy.sqrt((105.0 - c) / 140.0)
        t = numpy.sqrt(3.0 / 5.0)

        B1 = (77.0 - c) / 891.0
        B2 = (77.0 + c) / 891.0
        B3 = 25.0 / 324.0

        self.degree = 7
        data = [
            (B1, _symm_r_0(r)),
            (B2, _symm_r_0(s)),
            (B3, _pm2(t, t))
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0
        return
