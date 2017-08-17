# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, z, pm_array0, pm


class Ditkin(object):
    '''
    V. A. Ditkin,
    On certain approximate formulas for the calculation of triple integrals,
    Doklady Akad. Nauk SSSR (N.S.) 62 (1948), 445–447 (Russian).
    '''
    def __init__(self, index, alpha=0.0):
        if index == 1:
            self.degree = 5

            B0 = 4.0 / (alpha + 5.0)**2
            B1 = (alpha + 3.0) * (alpha + 7.0) / 12.0 / (alpha + 5.0)**2

            t = numpy.array([+1, -1])
            r, s = numpy.sqrt(
                    (alpha+5.0) * (5.0+t*numpy.sqrt(5.0)) / 10.0 / (alpha+7.0)
                    )

            data = [
                (B0, z(3)),
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                ]
        elif index == 2:
            self.degree = 5

            B0 = 4.0 / 25.0
            B1 = 21.0 / 500.0

            plus_minus = numpy.array([+1, -1])
            r, s = numpy.sqrt((15.0 + 5.0*plus_minus*numpy.sqrt(5.0)) / 42.0)
            t = numpy.sqrt(5.0 / 21.0)

            data = [
                (B0, z(3)),
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                (B1, pm(3, t)),
                ]
        else:
            assert index == 3, 'Illegal Ditkin index {}.'.format(index)
            self.degree = 7

            B0 = 16.0 / 175.0
            B1 = 81.0 / 1400.0
            B2 = 3.0 / 280.0

            plus_minus = numpy.array([+1, -1])
            sqrt5 = numpy.sqrt(5.0)
            r, s = numpy.sqrt((5.0 + plus_minus*sqrt5) / 18.0)
            t = numpy.sqrt(1.0/3.0)
            u, v = numpy.sqrt((3.0 - plus_minus*sqrt5) / 6.0)

            data = [
                (B0, z(3)),
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                (B2, pm_array0(3, [u, v], [0, 1])),
                (B2, pm_array0(3, [u, v], [1, 2])),
                (B2, pm_array0(3, [u, v], [2, 0])),
                (B2, pm(3, t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0/3.0 * numpy.pi
        return
