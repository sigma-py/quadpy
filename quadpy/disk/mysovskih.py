# -*- coding: utf-8 -*-
#
from math import sqrt, pi
import numpy

from ..helpers import untangle, z, fsd, fs_array


class Mysovskih(object):
    '''
    I.P. Mysovskih,
    On the construction of cubature formulas for the simplest regions,
    Z. Vychisl. Mat. i. Mat. Fiz. 4, 3-14, 1964.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index, alpha=0.0):
        if index == 1:
            self.degree = 4
            i = numpy.arange(5)
            b = sqrt((alpha + 4.0)/(alpha + 6.0))

            r = b * numpy.cos(0.4*i*pi)
            s = b * numpy.sin(0.4*i*pi)
            x = numpy.column_stack([r, s])

            B0 = 4.0 / (alpha + 4)**2
            B1 = (alpha + 2.0)*(alpha + 6.0) / 5.0 / (alpha + 4.0)**2

            data = [
                (B0, z(2)),
                (B1, x),
                ]
        elif index == 2:
            self.degree = 11

            sqrt10 = sqrt(10.0)
            sqrt601 = sqrt(601.0)

            B1 = (857.0*sqrt601 + 12707.0) / 20736.0 / sqrt601
            B2 = 125.0 / 3456.0
            B3 = (857.0*sqrt601 - 12707.0) / 20736.0 / sqrt601
            B4 = (340.0 + 25*sqrt10) / 10368.0
            B5 = (340.0 - 25*sqrt10) / 10368.0

            r1 = sqrt((31.0 - sqrt601) / 60.0)
            r2 = sqrt(3.0/5.0)
            r3 = sqrt((31.0 + sqrt601) / 60.0)
            r4 = sqrt((10.0 - sqrt10) / 20.0)
            r5 = sqrt((10.0 + sqrt10) / 20.0)

            s4 = sqrt((10.0 - sqrt10) / 60.0)
            s5 = sqrt((10.0 + sqrt10) / 60.0)

            data = [
                (B1, fsd(2, r1, 1)),
                (B2, fsd(2, r2, 1)),
                (B3, fsd(2, r3, 1)),
                (B4, fs_array([r4, s4])),
                (B5, fs_array([r5, s5])),
                ]
        else:
            assert index == 3
            # This is is the same as Rabinowitz-Richter
            self.degree = 15

            t = numpy.array([+1, -1])
            sqrt21 = sqrt(21.0)
            sqrt1401 = sqrt(1401.0)

            # ERR Stroud lists 5096 instead of 4998 here
            A1, A2 = (4998.0 + t * 343 * sqrt21) / 253125.0

            # ERR Stroud is missing the +- here
            B1, B3 = \
                (1055603.0 * sqrt1401 + t * 26076047.0) / 43.2e6 / sqrt1401
            B2 = 16807.0 / 8.0e5

            rho1, rho2 = numpy.sqrt((21.0 - t * sqrt21) / 28.0)
            sigma1, sigma3 = numpy.sqrt((69.0 - t * sqrt1401) / 112.0)
            sigma2 = sqrt(5.0/7.0)

            tau1 = 0.252863797091230
            tau2 = 0.577728928444823
            tau3 = 0.873836956644882
            tau4 = 0.989746802511491

            C1 = 0.398811120280412e-1
            C2 = 0.348550570365141e-1
            C3 = 0.210840370156484e-1
            C4 = 0.531979391979623e-2

            ka = numpy.arange(1, 9)
            xa = numpy.column_stack([
                numpy.cos((2*ka-1) / 8.0 * numpy.pi),
                numpy.sin((2*ka-1) / 8.0 * numpy.pi),
                ])

            kb = numpy.arange(1, 5)
            xb = numpy.column_stack([
                numpy.cos((2*kb-1) / 4.0 * numpy.pi),
                numpy.sin((2*kb-1) / 4.0 * numpy.pi),
                ])

            kc = numpy.arange(1, 5)
            xc = numpy.column_stack([
                numpy.cos(kc / 2.0 * numpy.pi),
                numpy.sin(kc / 2.0 * numpy.pi),
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
