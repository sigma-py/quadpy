# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import _symm_r_0, _z, _symm_s_t
from ..helpers import untangle


class Maxwell(object):
    '''
    J.C. Maxwell,
    On Approximate Multiple Integration between Limits by Summation.
    In W. Niven (Ed.), The Scientific Papers of James Clerk Maxwell,
    Cambridge Library Collection - Physical Sciences, pp. 604-611.
    Cambridge: Cambridge University Press.
    First published in 1890.
    <https://doi.org/10.1017/CBO9780511710377.061>.
    '''
    def __init__(self):
        self.name = 'Maxwell'
        self.degree = 7

        r = sqrt(fr(12, 35))
        s, t = [sqrt((93 + i*3*sqrt(186)) / 155) for i in [+1, -1]]

        data = [
            (fr(1, 81), _z()),
            (fr(49, 324), _symm_r_0(r)),
            # ERR typo in Stroud: 648 vs 649
            (fr(31, 648), _symm_s_t(s, t))
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
