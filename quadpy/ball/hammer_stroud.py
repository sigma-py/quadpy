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
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == '11-3':
            self.degree = 3
            data = [
                (1.0/6.0, fsd(3, sqrt(0.6), 1)),
                ]
        elif index == '12-3':
            self.degree = 5
            alpha = sqrt(3.0/7.0)
            data = [
                (1.0/15.0, z(3)),
                (7.0/90.0, fsd(3, alpha, 1)),
                (7.0/180.0, fsd(3, alpha, 2)),
                ]
        elif index in ['14-3a', '14-3b']:
            self.degree = 5

            t = 1 if index == '14-3a' else -1

            sqrt14 = sqrt(14.0)

            # ERR The article falsely gives 0.50824... instead of 0.050824...
            a1 = 1.0/125.0 * (9 + t * 2*sqrt14)
            c1 = (71.0 - t * 12 * sqrt14) / 1000.0

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

            t = 1 if index == '15-3a' else - 1

            sqrt30 = sqrt(30.0)
            nu2 = (45.0 - t * sqrt30)/57.0
            xi2 = (18.0 + t * sqrt30)/42.0
            eta2 = 7.0 / (27.0 + t * 2*sqrt30)

            # The extract expressions are from Stroud's book.
            a1 = 1.0 / 63.0 / nu2**3
            b1 = 1.0 / 630.0 / xi2**3
            c1 = 1.0 / 2520.0 / eta2**3
            a0 = 1.0 - 6*a1 - 12*b1 - 8*c1

            data = [
                (a0, z(3)),
                (a1, fsd(3, sqrt(nu2), 1)),
                (b1, fsd(3, sqrt(xi2), 2)),
                (c1, pm(3, sqrt(eta2))),
                ]
        self.points, self.weights = untangle(data)
        self.weights *= 4.0/3.0 * pi
        return
