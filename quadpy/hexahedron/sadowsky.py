# -*- coding: utf-8 -*-
#
from sympy import Rational as fr, sqrt

from .helpers import fs_r00, fs_rr0, fs_rrs
from ..helpers import untangle


class Sadowsky(object):
    '''
    Michael Sadowsky,
    A Formula for Approximate Computation of a Triple Integral,
    The American Mathematical Monthly,
    Vol. 47, No. 8 (Oct., 1940), pp. 539-543,
    <https://doi.org/10.2307/2303834>.
    '''
    def __init__(self):
        self.degree = 5
        data = [
            (fr(91, 450), fs_r00(1)),
            (fr(-20, 225), fs_rr0(1)),
            (fr(8, 225), fs_rrs(sqrt(fr(5, 8)), 1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
