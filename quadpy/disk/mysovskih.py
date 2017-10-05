# -*- coding: utf-8 -*-
#
import numpy
from sympy import sqrt, pi, sin, cos, Rational as fr

from ..helpers import untangle, z, fsd, fs_array


class Mysovskih(object):
    '''
    I.P. Mysovskih,
    On the construction of cubature formulas for the simplest regions,
    Z. Vychisl. Mat. i. Mat. Fiz. 4, 3-14, 1964.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index, alpha=0):
        if index == 1:
            self.degree = 4
            b = sqrt(fr(alpha + 4, alpha + 6))

            x = numpy.column_stack([
                [b * cos(2*i*pi/5) for i in range(5)],
                [b * sin(2*i*pi/5) for i in range(5)],
                ])

            B0 = fr(4, (alpha + 4)**2)
            B1 = fr((alpha + 2)*(alpha + 6), 5 * (alpha + 4)**2)

            data = [
                (B0, z(2)),
                (B1, x),
                ]
        elif index == 2:
            self.degree = 11

            sqrt10 = sqrt(10)
            sqrt601 = sqrt(601)

            B1, B3 = [
                (857*sqrt601 + t*12707) / 20736 / sqrt601
                for t in [+1, -1]
                ]
            B2 = fr(125, 3456)
            B4, B5 = [(340 + t * 25*sqrt10) / 10368 for t in [+1, -1]]

            r1, r3 = [sqrt((31 - t*sqrt601) / 60) for t in [+1, -1]]
            r2 = sqrt(fr(3, 5))
            r4, r5 = [sqrt((10 - t*sqrt10) / 20) for t in [+1, -1]]

            s4, s5 = [sqrt((10 - t*sqrt10) / 60) for t in [+1, -1]]

            data = [
                (B1, fsd(2, (r1, 1))),
                (B2, fsd(2, (r2, 1))),
                (B3, fsd(2, (r3, 1))),
                (B4, fs_array([r4, s4])),
                (B5, fs_array([r5, s5])),
                ]
        else:
            assert index == 3
            # This is is the same as Rabinowitz-Richter
            self.degree = 15

            sqrt21 = sqrt(21)
            sqrt1401 = sqrt(1401)

            # ERR Stroud lists 5096 instead of 4998 here
            A1, A2 = [(4998 + t * 343 * sqrt21) / 253125 for t in [+1, -1]]

            # ERR Stroud is missing the +- here
            B1, B3 = [
                (1055603 * sqrt1401 + t * 26076047) / 43200000 / sqrt1401
                for t in [+1, -1]
                ]

            B2 = fr(16807, 800000)

            rho1, rho2 = [sqrt((21 - t * sqrt21) / 28) for t in [+1, -1]]
            sigma1, sigma3 = [sqrt((69 - t*sqrt1401) / 112) for t in [+1, -1]]
            sigma2 = sqrt(fr(5, 7))

            tau1 = 0.252863797091230
            tau2 = 0.577728928444823
            tau3 = 0.873836956644882
            tau4 = 0.989746802511491

            C1 = 0.398811120280412e-1
            C2 = 0.348550570365141e-1
            C3 = 0.210840370156484e-1
            C4 = 0.531979391979623e-2

            xa = numpy.column_stack([
                [cos((2*k-1) * pi / 8) for k in range(1, 9)],
                [sin((2*k-1) * pi / 8) for k in range(1, 9)],
                ])

            xb = numpy.column_stack([
                [cos((2*k-1) * pi / 4) for k in range(1, 5)],
                [sin((2*k-1) * pi / 4) for k in range(1, 5)],
                ])

            xc = numpy.column_stack([
                [cos(k*pi / 2) for k in range(1, 5)],
                [sin(k*pi / 2) for k in range(1, 5)],
                ])

            data = [
                (A1, rho1 * xa),
                (A2, rho2 * xa),
                (B1, sigma1 * xb),
                (B2, sigma2 * xb),
                (B3, sigma3 * xb),
                (C1, tau1 * xc),
                (C2, tau2 * xc),
                (C3, tau3 * xc),
                (C4, tau4 * xc),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
