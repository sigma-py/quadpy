# -*- coding: utf-8 -*-
#
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

            # ERR Stroud falsely lists 7494893.0 instead of 7494892.0
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

            # The values are solutions of
            # 11025*x^3 - 19020*x^2 + 9370*x - 1212 = 0
            sigma2 = [
                0.2034846565025736047446226,
                0.5642857550696463424828808,
                0.9573996564549909371262381
                ]

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                 1432433.0 / 18849024.0,
                 1075.0 / 31104.0,
                 521.0 / 25920.0,
                ])
            B = numpy.linalg.solve(A, b)

            B0 = 2615.0 / 43632.0
            C = 16807.0 / 933120.0

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
                (B[0], numpy.sqrt(sigma2[0])*rs),
                (B[1], numpy.sqrt(sigma2[1])*rs),
                (B[2], numpy.sqrt(sigma2[2])*rs),
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

            # The values are solutions of
            # 4960228*x^4- 10267740*x^3 + 6746490*x^2 - 1476540*x + 70425 = 0
            sigma2 = [
                0.06530799472126887796787083,
                0.3071576960209604329098429,
                0.7363522053494293702057347,
                0.9611958210197321662727617,
                ]

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                 57719.0 / 675000.0,
                 9427.0 / 270000.0,
                 193.0 / 9000.0,
                 113.0 / 7200.0,
                ])
            ws = numpy.linalg.solve(A, b)

            data = [
                (ws[0], numpy.sqrt(sigma2[0]) * s),
                (ws[1], numpy.sqrt(sigma2[1]) * s),
                (ws[2], numpy.sqrt(sigma2[2]) * s),
                (ws[3], numpy.sqrt(sigma2[3]) * s),
                (wt1, tau1 * t),
                (wt2, tau2 * t),
                ]
        else:
            assert index == 8
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

            # The values are solutions of
            # 160901628*x^4 - 364759920*x^3 + 274856190*x^2 - 76570340*x
            # + 6054195 = 0
            sigma2 = [
                0.12953711368067995173112945,
                0.3813917897409179499603420,
                0.7814925125193852533809405,
                0.9745532848992097110929944,
                ]

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                121827491812.0 / 1802182496625.0,
                48541.0 / 1666980.0,
                977.0 / 55566.0,
                671.0 / 52920.0,
                ])
            ws = numpy.linalg.solve(A, b)

            data = [
                (m0, z(2)),
                (ws[0], numpy.sqrt(sigma2[0]) * s),
                (ws[1], numpy.sqrt(sigma2[1]) * s),
                (ws[2], numpy.sqrt(sigma2[2]) * s),
                (ws[3], numpy.sqrt(sigma2[3]) * s),
                (wt1, tau1 * t),
                (wt2, tau2 * t),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= numpy.pi
        return
