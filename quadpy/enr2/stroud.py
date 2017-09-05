# -*- coding: utf-8 -*-
#
from __future__ import division

from math import sqrt, pi

import numpy

from .stenger import Stenger
from .stroud1967a import Stroud1967a
from .stroud1967b import Stroud1967b
from .stroud_secrest import StroudSecrest

from ..helpers import untangle, fsd, pm, pm_array0


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index):
        self.name = 'Stroud_Enr2({})'.format(index)
        self.dim = n
        if index == '3-1':
            self.set_data(StroudSecrest(n, 'I'))
        elif index == '3-2':
            self.set_data(StroudSecrest(n, 'III'))
        elif index == '5-1a':
            self.set_data(Stroud1967a(n, 'a'))
        elif index == '5-1b':
            self.set_data(Stroud1967a(n, 'b'))
        elif index == '5-2':
            self.set_data(StroudSecrest(n, 'IV'))
        elif index == '5-3':
            self.degree = 5

            r = sqrt((n+2) / 4)
            s = sqrt((n+2) / 2 / (n-2))
            A = 4 / (n+2)**2
            B = (n-2)**2 / 2**n / (n+2)**2

            data = [
                (A, fsd(n, (r, 1))),
                (B, pm(n, s)),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= sqrt(pi)**n
        elif index == '5-4':
            # spherical product Lobatto
            self.degree = 5

            B0 = 2 / (n+2)
            data = [
                (B0, numpy.full((1, n), 0.0)),
                ]
            for k in range(1, n+1):
                rk = sqrt((k+2) / 2)
                s = sqrt(1/2)
                arr = [rk] + (n-k) * [s]
                idx = numpy.arange(k-1, n)
                data += [
                    (2**(k-n) / (k+1) / (k+2), pm_array0(n, arr, idx))
                    ]

            self.points, self.weights = untangle(data)
            self.weights *= sqrt(pi)**n
        elif index in ['5-5a', '5-5b']:
            self.degree = 5

            p_m = +1 if index == '5-5a' else -1

            r = sqrt((n + 2 + p_m * (n-1) * sqrt(2*(n+2))) / (2*n))
            s = sqrt((n + 2 - p_m * sqrt(2*(n+2))) / (2*n))
            A = 2 / (n+2)
            B = 1 / 2**n / (n+2)

            data = [
                (A, numpy.full((1, n), 0.0)),
                (B, fsd(n, (r, 1), (s, n-1))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= sqrt(pi)**n
        elif index == '5-6':
            assert n >= 5
            self.degree = 5

            sqrt2 = sqrt(2)
            sqrt2n1 = sqrt(2*(n+1))
            r = sqrt((n - sqrt2 + (n-1) * sqrt2n1) / (2*n))
            s = sqrt((n - sqrt2 - sqrt2n1) / (2*n))
            t = sqrt((1 + sqrt2) / 2)
            A = 1 / 2**n / (n+1)

            data = [
                (A, fsd(n, (r, 1), (s, n-1))),
                (A, pm(n, t)),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= sqrt(pi)**n
        elif index == '7-1a':
            self.set_data(Stroud1967b('2a', n))
        elif index == '7-1b':
            self.set_data(Stroud1967b('2b', n))
        elif index == '7-2':
            self.set_data(Stroud1967b('4', n))
        elif index == '7-3a':
            self.set_data(Stenger(n, 7, 'a'))
        elif index == '7-3b':
            self.set_data(Stenger(n, 7, 'b'))
        elif index == '9-1a':
            self.set_data(Stenger(n, 9, 'a'))
        elif index == '9-1b':
            self.set_data(Stenger(n, 9, 'b'))
        elif index == '11-1a':
            self.set_data(Stenger(n, 11, 'a'))
        else:
            assert index == '11-1b'
            self.set_data(Stenger(n, 11, 'b'))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
