# -*- coding: utf-8 -*-
#
"""
Shayne Waldron,
Symmetries of linear functionals,
University of Wisconsin-Madison, Center for Mathematical Sciences,
CMS Technical Summary Report #95-04, October 1994,
<http://ftp.cs.wisc.edu/Approx/symmetries.pdf>.
"""
from __future__ import division

import numpy

from .helpers import QuadrilateralScheme


def Waldron(r, alpha):
    assert r ** 2 >= 1 / 3

    R = r / numpy.sqrt(3 * r ** 2 - 1)

    beta = alpha + numpy.pi / 2

    points = numpy.array(
        [
            [+r * numpy.cos(alpha), +r * numpy.sin(alpha)],
            [-r * numpy.cos(alpha), -r * numpy.sin(alpha)],
            [+R * numpy.cos(beta), +R * numpy.sin(beta)],
            [-R * numpy.cos(beta), -R * numpy.sin(beta)],
        ]
    )

    weights = numpy.array(
        [2 / 3 / r ** 2, 2 / 3 / r ** 2, 2 / 3 / R ** 2, 2 / 3 / R ** 2]
    )
    return QuadrilateralScheme("Waldron", 3, weights, points)
