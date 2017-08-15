# -*- coding: utf-8 -*-
#
from math import pi, sqrt

from ..helpers import untangle, fsd, pm, z


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, index):
        if index == '11-3':
            self.degree = 3
            data = [
                (2.0/9.0, fsd(3, sqrt(0.6), 1)),
                ]
        elif index == '12-3':
            self.degree = 5
            alpha = sqrt(3.0/7.0)
            data = [
                (4.0/45.0, z(3)),
                (14.0/135.0, fsd(3, alpha, 1)),
                (7.0/135.0, fsd(3, alpha, 2)),
                ]
        elif index in ['14-3a', '14-3b']:
            self.degree = 5
            sqrt14 = sqrt(14.0)
            # TODO merge with +-
            if index == '14-3a':
                a1 = 4.0/375.0 * (9 + 2*sqrt14)
                c1 = (71.0 - 12 * sqrt14) / 750.0
                nu = sqrt((7.0 - sqrt14) / 7.0)
                eta1 = sqrt(5.0 / (21.0 - 2*sqrt14))
            else:
                assert index == '14-3b'
                a1 = 8.0/75.0 * (9.0 - 2*sqrt14)
                c1 = (71.0 + 12 * sqrt14) / 750.0
                nu = sqrt((7.0 + sqrt14) / 7.0)
                eta1 = sqrt(5.0 / (21.0 + 2*sqrt14))
            data = [
                (a1, fsd(3, nu, 1)),
                (c1, pm(3, eta1)),
                ]
        else:
            assert index in ['15-3a', '15-3b']
            # TODO continue here
            assert False

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
