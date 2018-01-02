# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy
import numpy

from ..helpers import (
    untangle, fsd, pm, z, pm_array, pm_array0, fs_array as fs
    )


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        pi = sympy.pi if symbolic else numpy.pi
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pm_ = numpy.array([+1, -1])

        if index == '11-2':
            self.degree = 3
            # ERR Wrongly stated in Stroud with 0.5 instead of sqrt(0.5)
            data = [
                (frac(1, 4), fsd(2, (sqrt(frac(1, 2)), 1))),
                ]
        elif index == '12-2':
            self.degree = 5
            data = [
                (frac(1, 6), z(2)),
                (frac(1, 6), fsd(2, (sqrt(frac(1, 2)), 1))),
                (frac(1, 24), pm(2, sqrt(frac(1, 2)))),
                ]
        elif index == '13-2':
            self.degree = 7
            sqrt29 = sqrt(29)
            b1, b2 = (551 + pm_ * 41 * sqrt29) / 6264
            xi1, xi2 = sqrt(3 / (9 + pm_*sqrt29) / 2)
            data = [
                (frac(2, 27), fsd(2, (sqrt(frac(3, 4)), 1))),
                (b1, pm(2, xi1)),
                (b2, pm(2, xi2)),
                ]
        elif index == '17':
            self.degree = 5
            data = [
                (frac(1, 4), z(2)),
                (frac(1, 8), pm_array0(2, [sqrt(frac(2, 3))], [0])),
                (frac(1, 8), pm_array([sqrt(frac(1, 6)), sqrt(frac(1, 2))])),
                ]
        elif index == '18':
            self.degree = 7
            data = [
                (frac(1, 16), fs([0.4247082002778669, 0.1759198966061612])),
                (frac(1, 16), fs([0.8204732385702833, 0.3398511429799874])),
                ]
        elif index == '19':
            self.degree = 9

            sqrt6 = sqrt(6)
            alpha1 = (16 + sqrt6) / 288
            alpha2 = (137 - 32*sqrt6) / 1818
            alpha3 = (520 + 155*sqrt6) / 3636 / 8

            data = [
                (frac(1, 9), z(2)),
                (alpha1, fs([0.5505043204538557, 0.2280263556769715])),
                (alpha2, fsd(2, (0.9192110607898046, 1))),
                (alpha3, fs([0.7932084745126058, 0.4645097310495256])),
                ]
        elif index == '20':
            self.degree = 11

            sqrt15 = sqrt(15)
            alpha1 = (34 - 5*sqrt15) / 396
            alpha2 = frac(5, 792) * (2 + sqrt15)
            alpha3 = 0.0727157433213629 / pi

            x1 = sqrt((5 - sqrt15) / 10)
            x2 = sqrt(frac(1, 2)) * x1
            x3 = sqrt((5 + sqrt15) / 10)

            data = [
                (frac(5, 144), fsd(2, (x1, 1))),
                (frac(5, 144), pm(2, x2)),
                (alpha1, fsd(2, (sqrt(0.5), 1))),
                (alpha2, fs([0.6125369400823741, 0.3532683074300921])),
                (alpha3, fs([0.8157480497746617, 0.4710132205252606])),
                (0.0727346698565653/pi, fsd(2, (x3, 1))),
                ]
        else:
            assert index == '21'
            self.degree = 15

            alpha0 = 0.0341505695624825 / pi
            alpha1 = 0.0640242008621985 / pi
            alpha2 = 0.0341505695624825 / pi

            data = [
                (alpha0, fs([0.2584361661674054, 0.0514061496288813])),
                (alpha1, fs([0.5634263397544869, 0.1120724670846205])),
                (alpha1, fs([0.4776497869993547, 0.3191553840796721])),
                (alpha1, fs([0.8028016728473508, 0.1596871812824163])),
                (alpha1, fs([0.6805823955716280, 0.4547506180649039])),
                (alpha2, fs([0.2190916025980981, 0.1463923286035535])),
                (alpha2, fs([0.9461239423417719, 0.1881957532057769])),
                (alpha2, fs([0.8020851487551318, 0.5359361621905023])),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
