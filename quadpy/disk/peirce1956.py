# -*- coding: utf-8 -*-
#
from __future__ import division
from math import sqrt, cos, sin, pi

from ..helpers import fsd, pm, untangle, fs_array


class Peirce1956(object):
    '''
    William Hollis Peirce,
    Numerical integration over planar regions,
    PhD thesis, University of Wisconsin--Madison, 1956,
    <https://books.google.de/books/about/Numerical_integration_over_planar_region.html?id=WR9SAAAAMAAJ&redir_esc=y>.

    '''
    def __init__(self, index):

        if index == 1:
            # Also: Formula 13-2 in Hammer-Stroud.
            self.degree = 7

            sqrt29 = sqrt(29.0)
            r = sqrt(0.75)
            s = sqrt((27 - 3 * sqrt29)/104.0)
            t = sqrt((27 + 3 * sqrt29)/104.0)

            B1 = 2.0/27.0
            # ERR Stroud falsely lists 4 instead of 41 here.
            B2 = (551.0 + 41*sqrt29) / 6264.0
            B3 = (551.0 - 41*sqrt29) / 6264.0

            data = [
                (B1, fsd(2, r, 1)),
                (B2, pm(2, s)),
                (B3, pm(2, t)),
                ]
        else:
            assert index == 2
            self.degree = 9

            sqrt15 = sqrt(15.0)
            cos_pi8 = cos(pi/8)
            sin_pi8 = sin(pi/8)

            r = sqrt((5 + sqrt15) / 10.0)
            u1 = sqrt((5.0 - sqrt15)/10.0) * cos_pi8
            v1 = sqrt((5.0 - sqrt15)/10.0) * sin_pi8
            u2 = sqrt(0.5) * cos_pi8
            v2 = sqrt(0.5) * sin_pi8
            u3 = sqrt((5 + sqrt15 + sqrt(185.0*sqrt15 - 700.0)) / 20.0)
            v3 = sqrt((5 + sqrt15 - sqrt(185.0*sqrt15 - 700.0)) / 20.0)

            B1 = (12060.0 - 1440.0*sqrt15) / 254088.0
            B2 = 5.0 / 144.0
            B3 = 1.0/18.0
            B4 = (5585.0 + 1440.0*sqrt15) / 508176.0

            data = [
                (B1, fsd(2, r, 1)),
                (B2, fs_array([u1, v1])),
                (B3, fs_array([u2, v2])),
                (B4, fs_array([u3, v3])),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
