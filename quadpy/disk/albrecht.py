# -*- coding: utf-8 -*-
#
import numpy
from mpmath import mp
from sympy import sin, cos, pi, Rational as fr, sqrt

from ..helpers import untangle, z, fsd, pm


class Albrecht(object):
    '''
    J. Albrecht,
    Formeln zur numerischen Integration über Kreisbereiche,
    Volume 40, Issue 10-11, 1960, Pages 514–517,
    <https://doi.org/10.1002/zamm.19600401014>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == 1:
            self.degree = 3

            t = numpy.column_stack([
                [cos((2*k+1)*pi / 4) for k in range(4)],
                [sin((2*k+1)*pi / 4) for k in range(4)],
                ])

            data = [
                (fr(1, 4), sqrt(fr(1, 2)) * t),
                ]
        elif index == 2:
            self.degree = 5

            t = numpy.array([
                [cos((2*k+1)*pi / 6) for k in range(6)],
                [sin((2*k+1)*pi / 6) for k in range(6)],
                ]).T

            data = [
                (fr(1, 4), z(2)),
                (fr(1, 8), sqrt(fr(2, 3)) * t),
                ]
        elif index == 3:
            self.degree = 7

            s = numpy.column_stack([
                [cos(2*k*pi / 4) for k in range(4)],
                [sin(2*k*pi / 4) for k in range(4)],
                ])
            t = numpy.column_stack([
                [cos((2*k+1)*pi / 4) for k in range(4)],
                [sin((2*k+1)*pi / 4) for k in range(4)],
                ])

            sqrt29 = sqrt(29)
            a1, a2 = [(551 + pm_ * 41 * sqrt29) / 6264 for pm_ in [+1, -1]]
            rho1, rho2 = [
                sqrt((27 - pm_ * 3 * sqrt29) / 52)
                for pm_ in [+1, -1]
                ]

            data = [
                (fr(2, 27), sqrt(fr(3, 4)) * t),
                (a1, rho1 * s),
                (a2, rho2 * s),
                ]
        elif index == 4:
            self.degree = 9

            sqrt111 = sqrt(111)
            rho1, rho2 = [
                sqrt((96 - pm_ * 4*sqrt111) / 155)
                for pm_ in [+1, -1]
                ]

            s = numpy.column_stack([
                [cos(2*k*pi / 6) for k in range(6)],
                [sin(2*k*pi / 6) for k in range(6)],
                ])

            t = numpy.column_stack([
                [cos((2*k+1)*pi / 6) for k in range(6)],
                [sin((2*k+1)*pi / 6) for k in range(6)],
                ])

            B0 = fr(251, 2304)
            B1, B2 = [
                (110297 + pm_ * 5713*sqrt111) / 2045952
                for pm_ in [+1, -1]
                ]
            C = fr(125, 3072)

            data = [
                (B0, z(2)),
                (B1, rho1 * s),
                (B2, rho2 * s),
                (C, sqrt(fr(4, 5)) * t),
                ]
        elif index == 5:
            self.degree = 11

            # The values are solutions of
            # 6317094x^3 - 10022245*x^2 + 4149900*x - 336375 = 0
            sigma2 = mp.polyroots([6317094, -10022245, 4149900, -336375])

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                fr(168899, 1350000),
                fr(7661, 180000),
                fr(71, 3000),
                ])
            B = mp.lu_solve(A, b)

            sqrt19 = sqrt(19)

            # ERR Stroud falsely lists sqrt(10) for s1.
            s1, s2 = [sqrt((125 - pm_ * 10*sqrt19) / 366) for pm_ in [+1, -1]]

            # ERR Stroud falsely lists 749489_3_.0 instead of 749489_2_.0
            C1, C2 = [
                (7494892 + pm_ * 1053263*sqrt19) / 205200000
                for pm_ in [+1, -1]
                ]
            D = fr(81, 3125)

            u = sqrt(fr(5, 6)) * cos(pi/8)
            v = sqrt(fr(5, 6)) * sin(pi/8)

            data = [
                (B[0], fsd(2, (sqrt(sigma2[0]), 1))),
                (B[1], fsd(2, (sqrt(sigma2[1]), 1))),
                (B[2], fsd(2, (sqrt(sigma2[2]), 1))),
                (C1, pm(2, s1)),
                (C2, pm(2, s2)),
                (D, fsd(2, (u, 1), (v, 1))),
                ]
        elif index == 6:
            self.degree = 13

            # The values are solutions of
            # 11025*x^3 - 19020*x^2 + 9370*x - 1212 = 0
            sigma2 = mp.polyroots([11025, -19020, 9370, -1212])

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                fr(1432433, 18849024),
                fr(1075, 31104),
                fr(521, 25920),
                ])
            B = mp.lu_solve(A, b)

            B0 = fr(2615, 43632)
            C = fr(16807, 933120)

            rs = numpy.column_stack([
                [cos(2*k*pi / 10) for k in range(10)],
                [sin(2*k*pi / 10) for k in range(10)],
                ])

            uv = numpy.column_stack([
                [cos((2*k+1)*pi / 10) for k in range(10)],
                [sin((2*k+1)*pi / 10) for k in range(10)],
                ])

            data = [
                (B0, z(2)),
                (B[0], sqrt(sigma2[0])*rs),
                (B[1], sqrt(sigma2[1])*rs),
                (B[2], sqrt(sigma2[2])*rs),
                (C, sqrt(fr(6, 7)) * uv)
                ]
        elif index == 7:
            self.degree = 15

            s = numpy.column_stack([
                [cos(2*k*pi / 8) for k in range(8)],
                [sin(2*k*pi / 8) for k in range(8)],
                ])

            t = numpy.column_stack([
                [cos((2*k+1)*pi / 8) for k in range(8)],
                [sin((2*k+1)*pi / 8) for k in range(8)],
                ])

            sqrt21 = sqrt(21)
            wt1, wt2 = [
                (4998 + pm_ * 343 * sqrt21) / 253125
                for pm_ in [+1, -1]
                ]
            tau1, tau2 = [sqrt((21 - pm_*sqrt21) / 28) for pm_ in [+1, -1]]

            # The values are solutions of
            # 4960228*x^4 - 10267740*x^3 + 6746490*x^2 - 1476540*x + 70425 = 0
            sigma2 = mp.polyroots([
                4960228, -10267740, 6746490, -1476540, 70425
                ])

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                fr(57719, 675000),
                fr(9427, 270000),
                fr(193, 9000),
                fr(113, 7200),
                ])
            ws = mp.lu_solve(A, b)

            data = [
                (ws[0], sqrt(sigma2[0]) * s),
                (ws[1], sqrt(sigma2[1]) * s),
                (ws[2], sqrt(sigma2[2]) * s),
                (ws[3], sqrt(sigma2[3]) * s),
                (wt1, tau1 * t),
                (wt2, tau2 * t),
                ]
        else:
            assert index == 8
            self.degree = 17

            s = numpy.column_stack([
                [cos(2*k*pi / 10) for k in range(10)],
                [sin(2*k*pi / 10) for k in range(10)],
                ])

            t = numpy.column_stack([
                [cos((2*k+1)*pi / 10) for k in range(10)],
                [sin((2*k+1)*pi / 10) for k in range(10)],
                ])

            m0 = fr(496439663, 13349499975)

            sqrt7 = sqrt(7)
            wt1, wt2 = [(125504 + pm_*16054*sqrt7)/8751645 for pm_ in [+1, -1]]
            tau1, tau2 = [sqrt((14 - pm_*sqrt7) / 18) for pm_ in [+1, -1]]

            # The values are solutions of
            # 160901628*x^4 - 364759920*x^3 + 274856190*x^2 - 76570340*x
            # + 6054195 = 0
            sigma2 = mp.polyroots([
                160901628, -364759920, 274856190, -76570340, 6054195
                ])

            A = numpy.vander(sigma2, increasing=True).T
            b = numpy.array([
                fr(121827491812, 1802182496625),
                fr(48541, 1666980),
                fr(977, 55566),
                fr(671, 52920),
                ])
            ws = mp.lu_solve(A, b)

            data = [
                (m0, z(2)),
                (ws[0], sqrt(sigma2[0]) * s),
                (ws[1], sqrt(sigma2[1]) * s),
                (ws[2], sqrt(sigma2[2]) * s),
                (ws[3], sqrt(sigma2[3]) * s),
                (wt1, tau1 * t),
                (wt2, tau2 * t),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
