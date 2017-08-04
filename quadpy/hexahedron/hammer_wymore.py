# -*- coding: utf-8 -*-
#
import numpy

from .helpers import fs_r00, fs_rr0, pm_rrr

from ..helpers import untangle


class HammerWymore(object):
    '''
    Preston C. Hammer and A. Wayne Wymore,
    Numerical evaluation of multiple integrals. I,
    Math. Comp. 11 (1957), 59-67,
    <https://doi.org/10.1090/S0025-5718-1957-0087220-6>.
    '''
    def __init__(self):
        self.degree = 7

        r2 = 6.0/7.0
        s2 = (960.0 - 3*numpy.sqrt(28798.0)) / 2726.0
        t2 = (960.0 + 3*numpy.sqrt(28798.0)) / 2726.0

        r = numpy.sqrt(r2)
        s = numpy.sqrt(s2)
        t = numpy.sqrt(t2)

        B1 = 1078.0 / 29160.0
        B2 = 343.0 / 29160.0
        B3 = (774.0*t2 - 230.0) / (9720.0 * (t2-s2))
        B4 = (230.0 - 774.*s2) / (9720.0 * (t2-s2))

        data = [
            (B1, fs_r00(r)),
            (B2, fs_rr0(r)),
            (B3, pm_rrr(s)),
            (B4, pm_rrr(t)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return
