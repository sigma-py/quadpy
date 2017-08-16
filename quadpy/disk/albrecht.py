# -*- coding: utf-8 -*-
#
import warnings

import numpy

from ..helpers import untangle, z, fsd, fsd2, pm


class Albrecht(object):
    '''
    J. Albrecht,
    Formeln zur numerischen Integration über Kreisbereiche,
    Volume 40, Issue 10-11, 1960, Pages 514–517,
    <https://dx.doi.org/10.1002/zamm.19600401014>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        t = numpy.array([1, -1])
        if index == 1:
            self.degree = 9

            sqrt111 = numpy.sqrt(111.0)
            rho1, rho2 = numpy.sqrt((96.0 - t * 4*sqrt111) / 155.0)

            k = numpy.arange(1, 7)
            x = numpy.array([
                numpy.cos(k*numpy.pi/3.0),
                numpy.sin(k*numpy.pi/3.0),
                ]).T

            uv = numpy.sqrt(0.8) * numpy.array([
                numpy.cos((2*k-1)*numpy.pi / 6.0),
                numpy.sin((2*k-1)*numpy.pi / 6.0),
                ]).T

            B0 = 251.0 / 2304.0
            B1, B2 = (110297.0 + t * 5713.0*sqrt111) / 2045952.0
            C = 125.0 / 3072.0

            data = [
                (B0, z(2)),
                (B1, rho1 * x),
                (B2, rho2 * x),
                (C, uv),
                ]
        else:
            assert index == 2
            warnings.warn('Albrecht\'s 2nd scheme is only single-precision.')
            self.degree = 11

            r1 = 0.326655862701
            r2 = 0.720984642976
            r3 = 0.979798373636

            B1 = 0.0655791415454
            B2 = 0.0490399916287
            B3 = 0.0104912371962

            sqrt19 = numpy.sqrt(19.0)

            # ERR Stroud falsely lists sqrt(10) for s1.
            s1, s2 = numpy.sqrt((125.0 - t * 10.0*sqrt19) / 366.0)

            C1, C2 = (7494893.0 + t * 1053263.0*sqrt19) / 205200000.0
            D = 81.0 / 3125.0

            u = numpy.sqrt(5.0/6.0) * numpy.cos(numpy.pi/8.0)
            v = numpy.sqrt(5.0/6.0) * numpy.sin(numpy.pi/8.0)

            data = [
                (B1, fsd(2, r1, 1)),
                (B2, fsd(2, r2, 1)),
                (B3, fsd(2, r3, 1)),
                (C1, pm(2, s1)),
                (C2, pm(2, s2)),
                (D, fsd2(2, u, v, 1, 1)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= numpy.pi
        return
