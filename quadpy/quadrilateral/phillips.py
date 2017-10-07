# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

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

        c = 3*sqrt(385)
        r, s = [sqrt((105 + i*c) / 140) for i in [+1, -1]]
        t = sqrt(fr(3, 5))

        B1, B2 = [(77 - i*c) / 891 for i in [+1, -1]]
        B3 = fr(25, 324)

        self.degree = 7
        data = [
            (B1, _symm_r_0(r)),
            (B2, _symm_r_0(s)),
            (B3, _pm2(t, t))
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
