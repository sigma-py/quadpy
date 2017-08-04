# -*- coding: utf-8 -*-
#
import numpy

from .helpers import fs_r00, fs_rr0, fs_rrs

from ..helpers import untangle


class Sadowsky(object):
    '''
    Michael Sadowsky,
    A Formula for Approximate Computation of a Triple Integral,
    The American Mathematical Monthly,
    Vol. 47, No. 8 (Oct., 1940), pp. 539-543,
    <https://dx.doi.org/10.2307/2303834>.
    '''
    def __init__(self):
        self.degree = 5
        data = [
            (91.0/450.0, fs_r00(1.0)),
            (-20.0/225.0, fs_rr0(1.0)),
            (8.0/225.0, fs_rrs(numpy.sqrt(5.0/8.0), 1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return
