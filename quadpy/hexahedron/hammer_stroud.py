# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, fsd, pm, z


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 5
            data = [
                (40.0/361.0, fsd(3, numpy.sqrt(19.0/30.0), 1)),
                (121.0/2888.0, pm(3, numpy.sqrt(19.0/33.0)))
                ]
        else:
            assert index in [2, 3]
            self.degree = 7

            i = 1.0 if index == 2 else -1.0

            r2 = (33.0 - i * numpy.sqrt(165.0)) / 28.0
            s2 = (30.0 + i * numpy.sqrt(165.0)) / 35.0
            t2 = (195.0 - i * 4.0*numpy.sqrt(165.0)) / 337.0

            r = numpy.sqrt(r2)
            s = numpy.sqrt(s2)
            t = numpy.sqrt(t2)

            B1 = 22.0/945.0 / r2**3
            B2 = 1.0/135.0 / s2**3
            B3 = 1.0/216.0 / t2**3
            B0 = 1.0 - 6.0*B1 - 12.0*B2 - 8.0*B3

            data = [
                (B0, z(3)),
                (B1, fsd(3, r, 1)),
                (B2, fsd(3, s, 2)),
                (B3, pm(3, t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return
