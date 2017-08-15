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

            t = 1 if index == '14-3a' else -1

            sqrt14 = sqrt(14.0)
            # ERR The article falsely gives 0.50824... instead of 0.050824...
            a1 = 4.0/375.0 * (9 + t * 2*sqrt14)
            c1 = (71.0 - t * 12 * sqrt14) / 750.0
            nu = sqrt((7.0 - t * sqrt14) / 7.0)
            eta1 = sqrt(5.0 / (21.0 - t * 2*sqrt14))

            data = [
                (a1, fsd(3, nu, 1)),
                (c1, pm(3, eta1)),
                ]
        else:
            assert index in ['15-3a', '15-3b'], \
                'Illegal index {}.'.format(index)

            self.degree = 7

            if index == '15-3a':
                a0 = 0.4156003482691997 / pi
                a1 = 0.1994483077968051 / pi
                b1 = 0.0380676101171267 / pi
                c1 = 0.2649610860413550 / pi
            else:
                assert index == '15-3b'
                a0 = 0.4441396821009518 / pi
                a1 = 0.0957384071760634 / pi
                b1 = 0.2508385364520637 / pi
                c1 = 0.0200197052755367 / pi

            t = 1 if index == '15-3a' else - 1

            sqrt30 = sqrt(30.0)
            nu = sqrt((45.0 - t * sqrt30)/57.0)
            xi = sqrt((18.0 + t * sqrt30)/42.0)
            eta = sqrt(7.0 / (27.0 + t * 2*sqrt30))

            data = [
                (a0, z(3)),
                (a1, fsd(3, nu, 1)),
                (b1, fsd(3, xi, 2)),
                (c1, pm(3, eta)),
                ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
