# -*- coding: utf-8 -*-
#
"""
G.M. Ewing,
On Approximate Cubature,
The American Mathematical Monthly,
Vol. 48, No. 2 (Feb., 1941), pp. 134-136,
<https://doi.org/10.2307/2303604>.
"""
from __future__ import division

import sympy

from .helpers import NCubeScheme
from ..helpers import untangle, z, pm


def Ewing(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [(frac(2, 3), z(n)), (frac(1, 3 * 2 ** n), pm(n, 1))]

    points, weights = untangle(data)
    reference_volume = 2 ** n
    weights *= reference_volume
    return NCubeScheme("Ewing", n, 3, weights, points)
