# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from ..helpers import untangle, fsd, pm, z


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, index):
        if index == '1-3':
            self.degree = 3
            data = [
                (fr(4, 3), fsd(3, (1, 1))),
                ]
        elif index == '2-3':
            self.degree = 5
            alpha = sqrt(fr(3, 5))
            data = [
                (+fr(56, 27), z(3)),
                (-fr(20, 81), fsd(3, (alpha, 1))),
                (+fr(50, 81), fsd(3, (alpha, 2))),
                ]
        elif index == '4-3':
            self.degree = 5
            data = [
                (fr(320, 361), fsd(3, (sqrt(fr(19, 30)), 1))),
                (fr(121, 361), pm(3, sqrt(fr(19, 33))))
                ]
        elif index in ['5-3a', '5-3b']:
            self.degree = 7

            i = 1 if index == '5-3a' else -1

            r2 = (33 - i * sqrt(165)) / 28
            s2 = (30 + i * sqrt(165)) / 35
            t2 = (195 - i * 4*sqrt(165)) / 337

            r = sqrt(r2)
            s = sqrt(s2)
            t = sqrt(t2)

            B1 = 176 / r2**3 / 945
            B2 = 8 / s2**3 / 135
            B3 = 8 / t2**3 / 216
            B0 = 8 - 6*B1 - 12*B2 - 8*B3

            data = [
                (B0, z(3)),
                (B1, fsd(3, (r, 1))),
                (B2, fsd(3, (s, 2))),
                (B3, pm(3, t)),
                ]
        else:
            assert index == '6-3'
            self.degree = 7
            alpha = sqrt(fr(6, 7))
            data = [
                (fr(1078, 3645), fsd(3, (alpha, 1))),
                (fr(343, 3645), fsd(3, (alpha, 2))),
                (0.2247031747656014, pm(3, 0.7341125287521153)),
                (0.4123338622714356, pm(3, 0.4067031864267161)),
                ]

        self.points, self.weights = untangle(data)
        return
