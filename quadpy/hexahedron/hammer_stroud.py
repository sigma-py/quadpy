# -*- coding: utf-8 -*-
#
import math

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
                (4.0/3.0, fsd(3, (1.0, 1))),
                ]
        elif index == '2-3':
            self.degree = 5
            alpha = math.sqrt(3.0/5.0)
            data = [
                (+56.0/27.0, z(3)),
                (-20.0/81.0, fsd(3, (alpha, 1))),
                (+50.0/81.0, fsd(3, (alpha, 2))),
                ]
        elif index == '4-3':
            self.degree = 5
            data = [
                (320.0/361.0, fsd(3, (math.sqrt(19.0/30.0), 1))),
                (121.0/361.0, pm(3, math.sqrt(19.0/33.0)))
                ]
        elif index in ['5-3a', '5-3b']:
            self.degree = 7

            i = 1.0 if index == '5-3a' else -1.0

            r2 = (33.0 - i * math.sqrt(165.0)) / 28.0
            s2 = (30.0 + i * math.sqrt(165.0)) / 35.0
            t2 = (195.0 - i * 4.0*math.sqrt(165.0)) / 337.0

            r = math.sqrt(r2)
            s = math.sqrt(s2)
            t = math.sqrt(t2)

            B1 = 176.0/945.0 / r2**3
            B2 = 8.0/135.0 / s2**3
            B3 = 8.0/216.0 / t2**3
            B0 = 8.0 - 6.0*B1 - 12.0*B2 - 8.0*B3

            data = [
                (B0, z(3)),
                (B1, fsd(3, (r, 1))),
                (B2, fsd(3, (s, 2))),
                (B3, pm(3, t)),
                ]
        else:
            assert index == '6-3'
            self.degree = 7
            alpha = math.sqrt(6.0/7.0)
            data = [
                (1078.0/3645.0, fsd(3, (alpha, 1))),
                (343.0/3645.0, fsd(3, (alpha, 2))),
                (0.2247031747656014, pm(3, 0.7341125287521153)),
                (0.4123338622714356, pm(3, 0.4067031864267161)),
                ]

        self.points, self.weights = untangle(data)
        return
