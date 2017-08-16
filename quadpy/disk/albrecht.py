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
        plus_minus = numpy.array([+1, -1])
        if index == 1:
            self.degree = 3

            k = numpy.arange(4)
            t = numpy.array([
                numpy.cos((2*k+1)*numpy.pi / 4.0),
                numpy.sin((2*k+1)*numpy.pi / 4.0),
                ]).T

            data = [
                (0.25, numpy.sqrt(0.5) * t),
                ]
        elif index == 2:
            self.degree = 5

            k = numpy.arange(6)
            t = numpy.array([
                numpy.cos((2*k+1)*numpy.pi / 6.0),
                numpy.sin((2*k+1)*numpy.pi / 6.0),
                ]).T

            data = [
                (0.25, z(2)),
                (0.125, numpy.sqrt(2.0/3.0) * t),
                ]
        elif index == 3:
            self.degree = 7

            k = numpy.arange(4)
            s = numpy.array([
                numpy.cos(2*numpy.pi * k/4.0),
                numpy.sin(2*numpy.pi * k/4.0),
                ]).T
            t = numpy.array([
                numpy.cos((2*k+1)*numpy.pi / 4.0),
                numpy.sin((2*k+1)*numpy.pi / 4.0),
                ]).T

            sqrt29 = numpy.sqrt(29.0)
            a1, a2 = (551.0 + plus_minus * 41 * sqrt29) / 6264.0
            rho1, rho2 = numpy.sqrt((27.0 - plus_minus * 3 * sqrt29) / 52.0)

            data = [
                (2.0/27.0, numpy.sqrt(3.0/4.0) * t),
                (a1, rho1 * s),
                (a2, rho2 * s),
                ]
        elif index == 4:
            self.degree = 9

            sqrt111 = numpy.sqrt(111.0)
            rho1, rho2 = numpy.sqrt((96.0 - plus_minus * 4*sqrt111) / 155.0)

            k = numpy.arange(6)
            s = numpy.array([
                numpy.cos(2*numpy.pi * k/6.0),
                numpy.sin(2*numpy.pi * k/6.0),
                ]).T

            t = numpy.array([
                numpy.cos((2*k+1)*numpy.pi / 6.0),
                numpy.sin((2*k+1)*numpy.pi / 6.0),
                ]).T

            B0 = 251.0 / 2304.0
            B1, B2 = (110297.0 + plus_minus * 5713.0*sqrt111) / 2045952.0
            C = 125.0 / 3072.0

            data = [
                (B0, z(2)),
                (B1, rho1 * s),
                (B2, rho2 * s),
                (C, numpy.sqrt(0.8) * t),
                ]
        elif index == 5:
            warnings.warn('Albrecht\'s scheme no. 5 is only single-precision.')
            self.degree = 11

            r1 = 0.326655862701
            r2 = 0.720984642976
            r3 = 0.979798373636

            B1 = 0.0655791415454
            B2 = 0.0490399916287
            B3 = 0.0104912371962

            sqrt19 = numpy.sqrt(19.0)

            # ERR Stroud falsely lists sqrt(10) for s1.
            s1, s2 = numpy.sqrt((125.0 - plus_minus * 10.0*sqrt19) / 366.0)

            C1, C2 = (7494893.0 + plus_minus * 1053263.0*sqrt19) / 205200000.0
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
        else:
            assert index == 6
            self.degree = 13

            B0 = 2615.0 / 43632.0
            B1 = 0.0314864413570
            B2 = 0.0367783672793
            B3 = 0.00773026675860
            C = 16807.0 / 933120.0

            rho1 = 0.451092736034
            rho2 = 0.751189560011
            rho3 = 0.978468015039

            k = numpy.arange(10)
            rs = numpy.array([
                numpy.cos(2*k*numpy.pi/10.0),
                numpy.sin(2*k*numpy.pi/10.0),
                ]).T

            uv = numpy.array([
                numpy.cos((2*k+1)*numpy.pi/10.0),
                numpy.sin((2*k+1)*numpy.pi/10.0),
                ]).T

            data = [
                (B0, z(2)),
                (B1, rho1*rs),
                (B2, rho2*rs),
                (B3, rho3*rs),
                (C, numpy.sqrt(6.0/7.0) * uv)
                ]

        self.points, self.weights = untangle(data)
        self.weights *= numpy.pi
        return
