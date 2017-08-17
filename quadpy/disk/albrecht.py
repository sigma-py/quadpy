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
            self.degree = 11

            # The values are solutions of
            # 6317094x^3- 10022245*x^2 + 4149900*x - 336375 = 0
            sigma2 = [
                0.10670405263689525737465523,
                0.5198188554069141664267916,
                0.9600048529804843416546315,
                ]

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                168899.0 / 1350000.0,
                7661.0 / 180000.0,
                71.0 / 3000.0,
                ])
            B = numpy.linalg.solve(A, b)

            sqrt19 = numpy.sqrt(19.0)

            # ERR Stroud falsely lists sqrt(10) for s1.
            s1, s2 = numpy.sqrt((125.0 - plus_minus * 10.0*sqrt19) / 366.0)

            # TODO check if stroud has 7494893.0 instead of 7494892.0
            C1, C2 = (7494892.0 + plus_minus * 1053263.0*sqrt19) / 205200000.0
            D = 81.0 / 3125.0

            u = numpy.sqrt(5.0/6.0) * numpy.cos(numpy.pi/8.0)
            v = numpy.sqrt(5.0/6.0) * numpy.sin(numpy.pi/8.0)

            data = [
                (B[0], fsd(2, numpy.sqrt(sigma2[0]), 1)),
                (B[1], fsd(2, numpy.sqrt(sigma2[1]), 1)),
                (B[2], fsd(2, numpy.sqrt(sigma2[2]), 1)),
                (C1, pm(2, s1)),
                (C2, pm(2, s2)),
                (D, fsd2(2, u, v, 1, 1)),
                ]
        elif index == 6:
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
        elif index == 7:
            self.degree = 15

            k = numpy.arange(8)
            s = numpy.array([
                numpy.cos(2*numpy.pi * k/8.0),
                numpy.sin(2*numpy.pi * k/8.0),
                ]).T

            t = numpy.array([
                numpy.cos((2*k+1)*numpy.pi / 8.0),
                numpy.sin((2*k+1)*numpy.pi / 8.0),
                ]).T

            sqrt21 = numpy.sqrt(21.0)
            wt1, wt2 = (4998 + plus_minus * 343 * sqrt21) / 253125.0
            tau1, tau2 = numpy.sqrt((21.0 - plus_minus * sqrt21) / 28.0)

            ws1 = 0.204136860290e-1
            ws2 = 0.371360833569e-1
            ws3 = 0.209029582465e-1
            ws4 = 0.705690199725e-2

            sigma1 = 0.255554289186
            sigma2 = 0.554218094274
            sigma3 = 0.858109669768
            sigma4 = 0.980405947054

            data = [
                (ws1, sigma1 * s),
                (ws2, sigma2 * s),
                (ws3, sigma3 * s),
                (ws4, sigma4 * s),
                (wt1, tau1 * t),
                (wt2, tau2 * t),
                ]
        else:
            assert index == 8
            warnings.warn('Albrecht\'s scheme no. 8 is only single-precision.')
            self.degree = 17

            k = numpy.arange(10)
            s = numpy.array([
                numpy.cos(2*numpy.pi * k/10.0),
                numpy.sin(2*numpy.pi * k/10.0),
                ]).T

            t = numpy.array([
                numpy.cos((2*k+1)*numpy.pi / 10.0),
                numpy.sin((2*k+1)*numpy.pi / 10.0),
                ]).T

            m0 = 496439663.0 / 13349499975.0

            sqrt7 = numpy.sqrt(7.0)
            wt1, wt2 = (125504.0 + plus_minus * 16054 * sqrt7) / 8751645.0
            tau1, tau2 = numpy.sqrt((14.0 - plus_minus * sqrt7) / 18.0)

            ws1 = 0.206024726860e-1
            ws2 = 0.277365659974e-1
            ws3 = 0.150158249601e-1
            ws4 = 0.424511227320e-2

            sigma1 = 0.359912647292
            sigma2 = 0.617569259064
            sigma3 = 0.884020651636
            sigma4 = 0.987194654007

            data = [
                (m0, z(2)),
                (ws1, sigma1 * s),
                (ws2, sigma2 * s),
                (ws3, sigma3 * s),
                (ws4, sigma4 * s),
                (wt1, tau1 * t),
                (wt2, tau2 * t),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= numpy.pi
        return
