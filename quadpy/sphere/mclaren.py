# -*- coding: utf-8 -*-
#
import math
import numpy

from ..helpers import untangle, pm_array0, fsd, pm_array, pm, fsd2


class McLaren(object):
    '''
    A.D. McLaren,
    Optimal Numerical Integration on a Sphere,
    Mathematics of Computation, Vol. 17, No. 84. (Oct., 1963), pp. 361-383,
    <https://doi.org/10.1090/S0025-5718-1963-0159418-2>.
    '''
    # <https://github.com/PyCQA/pylint/issues/1472>
    # pylint: disable=too-many-locals, invalid-unary-operand-type
    def __init__(self, index):
        if index == 1:
            self.degree = 3

            data = [
                (1.0/12.0, fsd(3, math.sqrt(0.5), 2))
                ]
        elif index == 2:
            self.degree = 5

            # Stroud doesn't mention u=1, but it's implied. (After all, this is
            # integration on a sphere.)
            u = 1.0

            r = 0.5
            plus_minus = numpy.array([+1, -1])
            s, t = (math.sqrt(5.0) + plus_minus * 1) / 4.0

            data = [
                (1.0/30.0, fsd(3, u, 1)),
                (1.0/30.0, pm_array([r, s, t])),
                (1.0/30.0, pm_array([t, r, s])),
                (1.0/30.0, pm_array([s, t, r])),
                ]
        elif index == 3:
            self.degree = 7

            # the positive roots of
            #  z^6 - z^4 + 0.2*z^2 - 1/105 = 0,
            # i.e., the squares of the roots of
            #  z^3 - z^2 + 0.2*z^1 - 1/105 = 0,
            r = math.sqrt(0.7503835498819236124134653)
            s = math.sqrt(0.1785220127761020457111678)
            t = math.sqrt(0.07109443734197434187536689)

            u = numpy.array([+r, -r, +s, -s, +t, -t])
            v = numpy.array([+s, +t, +t, +r, +r, +s])
            w = numpy.array([+t, +s, +r, +t, +s, +r])

            data = [
                (1.0/24.0, numpy.column_stack([+u, +v, +w])),
                (1.0/24.0, numpy.column_stack([+u, -v, -w])),
                (1.0/24.0, numpy.column_stack([+u, +w, -v])),
                (1.0/24.0, numpy.column_stack([+u, -w, +v])),
                ]
        elif index == 4:
            self.degree = 8

            # the positive roots of
            #  z^6 - z^4 + 5/21 * z^2 - 5/441 = 0,
            # i.e., the squares of the roots of
            #  z^3 - z^2 + 5/21 * z^1 - 5/441 = 0,
            r = math.sqrt(0.6697999083949878692486451)
            s = math.sqrt(0.2667404813476434816325695)
            t = math.sqrt(0.06345961025736864911878537)

            u = numpy.array([+r, -r, +s, -s, +t, -t])
            v = numpy.array([+s, +t, +t, +r, +r, +s])
            w = numpy.array([+t, +s, +r, +t, +s, +r])

            data = [
                (16.0/600.0, fsd(3, 1.0, 1)),
                (21.0/600.0, numpy.column_stack([+u, +v, +w])),
                (21.0/600.0, numpy.column_stack([+u, -v, -w])),
                (21.0/600.0, numpy.column_stack([+u, +w, -v])),
                (21.0/600.0, numpy.column_stack([+u, -w, +v])),
                ]
        elif index == 5:
            self.degree = 9

            plus_minus = numpy.array([+1, -1])
            r, s = numpy.sqrt((5.0 + plus_minus * math.sqrt(5.0)) / 10.0)
            u, v = numpy.sqrt((3.0 - plus_minus * math.sqrt(5.0)) / 6.0)
            t = math.sqrt(1.0/3.0)

            B1 = 25.0 / 840.0
            B2 = 27.0 / 840.0

            data = [
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                #
                (B2, pm_array0(3, [u, v], [0, 1])),
                (B2, pm_array0(3, [u, v], [1, 2])),
                (B2, pm_array0(3, [u, v], [2, 0])),
                #
                (B2, pm(3, t)),
                ]
        elif index == 6:
            self.degree = 9

            plus_minus = numpy.array([+1, -1])
            r, s = numpy.sqrt((5.0 + plus_minus * math.sqrt(5.0)) / 10.0)
            t = 1.0
            u = 0.5
            v, w = (math.sqrt(5.0) + plus_minus * 1) / 4.0

            B = 25.0 / 1260.0
            C = 32.0 / 1260.0

            data = [
                # ERR Stroud is missing +- at the first r.
                (B, pm_array0(3, [r, s], [0, 1])),
                (B, pm_array0(3, [r, s], [1, 2])),
                (B, pm_array0(3, [r, s], [2, 0])),
                #
                (C, fsd(3, t, 1)),
                #
                (C, pm_array([u, v, w])),
                (C, pm_array([w, u, v])),
                (C, pm_array([v, w, u])),
                ]
        elif index == 7:
            self.degree = 9

            plus_minus = numpy.array([+1, -1])
            r, s = numpy.sqrt((3.0 - plus_minus * math.sqrt(5.0)) / 6.0)
            t = math.sqrt(1.0/3.0)
            # ERR Stroud falsely gives sqrt(0.5)
            u = 0.5
            v, w = (math.sqrt(5.0) + plus_minus * 1) / 4.0

            B = -9.0 / 140.0
            C = 16.0 / 210.0

            data = [
                (B, pm_array0(3, [r, s], [0, 1])),
                (B, pm_array0(3, [r, s], [1, 2])),
                (B, pm_array0(3, [r, s], [2, 0])),
                #
                (B, pm(3, t)),
                #
                (C, fsd(3, 1.0, 1)),
                #
                (C, pm_array([u, v, w])),
                (C, pm_array([w, u, v])),
                (C, pm_array([v, w, u])),
                ]
        elif index == 8:
            self.degree = 11

            r = 1.0
            s = math.sqrt(0.5)
            t = math.sqrt(1.0/3.0)

            u = math.sqrt(1.0/11.0)
            v = math.sqrt(9.0/11.0)

            B1 = 9216.0 / 725760.0
            B2 = 16384.0 / 725760.0
            B3 = 15309.0 / 725760.0
            B4 = 14641.0 / 725760.0

            data = [
                (B1, fsd(3, r, 1)),
                (B2, fsd(3, s, 2)),
                (B3, pm(3, t)),
                (B4, fsd2(3, u, v, 2, 1)),
                ]
        elif index == 9:
            self.degree = 11

            plus_minus = numpy.array([+1, -1])
            sqrt5 = math.sqrt(5.0)

            p, q = numpy.sqrt((5.0 + plus_minus*sqrt5) / 10.0)
            r, s = numpy.sqrt((3.0 - plus_minus*sqrt5) / 6.0)
            t = math.sqrt(1.0/3.0)

            u = 0.5
            v, w = (math.sqrt(5) + plus_minus * 1) / 4.0

            B = 625.0 / 27720.0
            C = 243.0 / 27720.0
            D = 512.0 / 27720.0

            data = [
                (B, pm_array0(3, [p, q], [0, 1])),
                (B, pm_array0(3, [p, q], [1, 2])),
                (B, pm_array0(3, [p, q], [2, 0])),
                #
                (C, pm_array0(3, [r, s], [0, 1])),
                (C, pm_array0(3, [r, s], [1, 2])),
                (C, pm_array0(3, [r, s], [2, 0])),
                #
                (C, pm(3, t)),
                #
                (D, fsd(3, 1.0, 1)),
                #
                (D, pm_array([u, v, w])),
                (D, pm_array([w, u, v])),
                (D, pm_array([v, w, u])),
                ]
        else:
            assert index == 10
            self.degree = 14

            plus_minus = numpy.array([+1, -1])

            r, s = numpy.sqrt((5.0 - plus_minus * numpy.sqrt(5.0)) / 10.0)
            B = 125.0 / 10080.0
            C = 143.0 / 10080.0

            # The roots of
            #
            # 2556125 y^6 - 5112250 y^5 + 3578575 y^4 - 1043900 y^3
            #     + 115115 y^2 - 3562 y + 9 =0
            #
            # in decreasing order.
            y = [
                0.8318603575087328951583062165711519728388,
                0.5607526046766541293084396308069013490725,
                0.4118893592345073860321480490176804941547,
                0.1479981814629634692260834719469411619893,
                0.04473134613410273910111648293922113227845,
                0.002768150983039381173906148718103889666260,
                ]
            z = numpy.sqrt(y)

            u = numpy.array([
                z[3] - z[2],
                z[1] - z[4],
                z[5] - z[1],
                z[2] - z[5],
                z[4] - z[3],
                ]) / 2.0 / s
            v = numpy.array([
                z[4] + z[5],
                z[5] + z[3],
                z[2] + z[4],
                z[3] + z[1],
                z[1] + z[2],
                ]) / 2.0 / s
            w = numpy.array([
                z[0] + z[1],
                z[0] + z[2],
                z[0] + z[3],
                z[0] + z[4],
                z[0] + z[5],
                ]) / 2.0 / s

            data = [
                (B, pm_array0(3, [r, s], [0, 1])),
                (B, pm_array0(3, [r, s], [1, 2])),
                (B, pm_array0(3, [r, s], [2, 0])),
                #
                (C, numpy.column_stack([+u, +v, +w])),
                (C, numpy.column_stack([+u, -v, -w])),
                (C, numpy.column_stack([-u, -v, +w])),
                (C, numpy.column_stack([-u, +v, -w])),
                #
                (C, numpy.column_stack([+v, +w, +u])),
                (C, numpy.column_stack([+v, -w, -u])),
                (C, numpy.column_stack([-v, -w, +u])),
                (C, numpy.column_stack([-v, +w, -u])),
                #
                (C, numpy.column_stack([+w, +u, +v])),
                (C, numpy.column_stack([+w, -u, -v])),
                (C, numpy.column_stack([-w, -u, +v])),
                (C, numpy.column_stack([-w, +u, -v])),
                ]

        self.points, self.weights = untangle(data)
        return
