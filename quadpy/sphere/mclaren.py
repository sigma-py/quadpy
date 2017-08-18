# -*- coding: utf-8 -*-
#
import math
import numpy

from ..helpers import untangle, pm_array0, fsd


class McLaren(object):
    '''
    A.D. McLaren,
    Optimal Numerical Integration on a Sphere,
    Mathematics of Computation, Vol. 17, No. 84. (Oct., 1963), pp. 361-383,
    <https://doi.org/10.1090/S0025-5718-1963-0159418-2>.
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 3

            data = [
                (1.0/12.0, fsd(3, math.sqrt(0.5), 2))
                ]
        else:
            assert index == 9
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
