# -*- coding: utf-8 -*-
#
from math import pi, sqrt

from ..helpers import untangle, fsd, pm, z, pm_array, pm_array0, fs_array


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == '11-2':
            self.degree = 3
            # ERR Wrongly stated in Stroud with 0.5 instead of sqrt(0.5)
            data = [
                (0.25, fsd(2, (sqrt(0.5), 1))),
                ]
        elif index == '12-2':
            self.degree = 5
            data = [
                (1.0/6.0, z(2)),
                (1.0/6.0, fsd(2, (sqrt(0.5), 1))),
                (1.0/24.0, pm(2, sqrt(0.5))),
                ]
        elif index == '13-2':
            self.degree = 7
            sqrt29 = sqrt(29.0)
            b1 = (551.0 + 41.0 * sqrt29) / 6264.0
            b2 = (551.0 - 41.0 * sqrt29) / 6264.0
            xi1 = sqrt(3.0 / 2 / (9 + sqrt29))
            xi2 = sqrt(3.0 / 2 / (9 - sqrt29))
            data = [
                (2.0/27.0, fsd(2, (sqrt(0.75), 1))),
                (b1, pm(2, xi1)),
                (b2, pm(2, xi2)),
                ]
        elif index == '17':
            self.degree = 5
            data = [
                (1.0/4.0, z(2)),
                (1.0/8.0, pm_array0(2, [sqrt(2.0/3.0)], [0])),
                (1.0/8.0, pm_array([sqrt(1.0/6.0), sqrt(1.0/2.0)])),
                ]
        elif index == '18':
            self.degree = 7
            data = [
                (1.0/16.0, fs_array([0.4247082002778669, 0.1759198966061612])),
                (1.0/16.0, fs_array([0.8204732385702833, 0.3398511429799874])),
                ]
        elif index == '19':
            self.degree = 9

            sqrt6 = sqrt(6.0)
            alpha1 = (16.0 + sqrt6) / 288.0
            alpha2 = (137.0 - 32*sqrt6) / 1818.0
            alpha3 = (520.0 + 155*sqrt6) / 3636.0 / 8

            data = [
                (1.0/9.0, z(2)),
                (alpha1, fs_array([0.5505043204538557, 0.2280263556769715])),
                (alpha2, fsd(2, (0.9192110607898046, 1))),
                (alpha3, fs_array([0.7932084745126058, 0.4645097310495256])),
                ]
        elif index == '20':
            self.degree = 11

            sqrt15 = sqrt(15.0)
            alpha1 = (34.0 - 5*sqrt15) / 396.0
            alpha2 = 5.0/792.0 * (2 + sqrt15)
            alpha3 = 0.0727157433213629 / pi

            x1 = sqrt((5.0 - sqrt15) / 10.0)
            x2 = sqrt(0.5) * x1
            x3 = sqrt((5.0 + sqrt15) / 10.0)

            data = [
                (5.0/144.0, fsd(2, (x1, 1))),
                (5.0/144.0, pm(2, x2)),
                (alpha1, fsd(2, (sqrt(0.5), 1))),
                (alpha2, fs_array([0.6125369400823741, 0.3532683074300921])),
                (alpha3, fs_array([0.8157480497746617, 0.4710132205252606])),
                (0.0727346698565653/pi, fsd(2, (x3, 1))),
                ]
        else:
            assert index == '21'
            self.degree = 15

            alpha0 = 0.0341505695624825 / pi
            alpha1 = 0.0640242008621985 / pi
            alpha2 = 0.0341505695624825 / pi

            data = [
                (alpha0, fs_array([0.2584361661674054, 0.0514061496288813])),
                (alpha1, fs_array([0.5634263397544869, 0.1120724670846205])),
                (alpha1, fs_array([0.4776497869993547, 0.3191553840796721])),
                (alpha1, fs_array([0.8028016728473508, 0.1596871812824163])),
                (alpha1, fs_array([0.6805823955716280, 0.4547506180649039])),
                (alpha2, fs_array([0.2190916025980981, 0.1463923286035535])),
                (alpha2, fs_array([0.9461239423417719, 0.1881957532057769])),
                (alpha2, fs_array([0.8020851487551318, 0.5359361621905023])),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
