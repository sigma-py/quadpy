# -*- coding: utf-8 -*-
#
import numpy

from .helpers import fs_r00, pm_rrr, fs_rr0, z

from ..helpers import untangle


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
                (40.0/361.0, fs_r00(numpy.sqrt(19.0/30.0))),
                (121.0/2888.0, pm_rrr(numpy.sqrt(19.0/33.0)))
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
                (B0, z()),
                (B1, fs_r00(r)),
                (B2, fs_rr0(s)),
                (B3, pm_rrr(t))
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return
