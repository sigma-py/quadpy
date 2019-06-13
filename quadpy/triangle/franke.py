# -*- coding: utf-8 -*-
#
"""
Richard Franke,
Obtaining cubatures for rectangles and other planar regions by using orthogonal
polynomials,
Math. Comp. 25 (1971), 803-817,
<https://doi.org/10.1090/S0025-5718-1971-0300440-5>.
"""
from __future__ import division

import numpy

from .helpers import TriangleScheme


def franke_9():
    a1 = 0.646341098016171e-1
    a2 = 0.250478764260821
    a3 = 0.405288113134598
    a4 = 0.483428507060240
    b1 = 0.490241549057468e-1
    c1 = 0.312418129002285
    b2 = 0.272654917225016e-1
    c2 = 0.649829918830148
    b3 = 0.748092005042521e-2
    c3 = 0.922929224698637
    b4 = 0.166718687651425
    c4 = 0.775796880494268
    b5 = 0.151969575382297
    c5 = 0.569101341800312

    points = numpy.array(
        [
            [a1, a1],
            [a2, a2],
            [a3, a3],
            [a4, a4],
            [b1, c1],
            [c1, b1],
            [b2, c2],
            [c2, b2],
            [b3, c3],
            [c3, b3],
            [b4, c4],
            [c4, b4],
            [b5, c5],
            [c5, b5],
        ]
    )
    weights = 2 * numpy.array(
        [
            0.263321501360460e-1,
            0.666750609902085e-1,
            0.598398472297514e-1,
            0.302244308027287e-1,
            0.387139102462897e-1,
            0.387139102462897e-1,
            0.223103130816147e-1,
            0.223103130816147e-1,
            0.930956404694027e-2,
            0.930956404694027e-2,
            0.365382927009296e-1,
            0.365382927009296e-1,
            0.515921753448585e-1,
            0.515921753448585e-1,
        ]
    )
    bary = numpy.array([points[:, 0], points[:, 1], 1 - numpy.sum(points, axis=1)]).T

    return TriangleScheme("Franke 8", 7, weights, points, bary)


Franke = {"9": franke_9}
