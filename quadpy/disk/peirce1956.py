# -*- coding: utf-8 -*-
#
from __future__ import division
import math

from ..helpers import fsd, pm, untangle


class Peirce1956(object):
    '''
    William Hollis Peirce,
    Numerical integration over planar regions,
    PhD thesis, University of Wisconsin--Madison, 1956,
    <https://books.google.de/books/about/Numerical_integration_over_planar_region.html?id=WR9SAAAAMAAJ&redir_esc=y>.

    Also: Formula 13-2 in Hammer-Stroud.
    '''
    def __init__(self):
        self.degree = 7

        sqrt29 = math.sqrt(29.0)
        r = math.sqrt(0.75)
        s = math.sqrt((27 - 3 * sqrt29)/104.0)
        t = math.sqrt((27 + 3 * sqrt29)/104.0)

        B1 = 2.0/27.0
        # ERR Stroud falsely lists 4 instead of 41 here.
        B2 = (551.0 + 41*sqrt29) / 6264.0
        B3 = (551.0 - 41*sqrt29) / 6264.0

        data = [
            (B1, fsd(2, r, 1)),
            (B2, pm(2, s)),
            (B3, pm(2, t)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= math.pi
        return
