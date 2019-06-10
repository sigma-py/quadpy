# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy


class Waldron(object):
    """
    Shayne Waldron,
    Symmetries of linear functionals,
    University of Wisconsin-Madison, Center for Mathematical Sciences,
    CMS Technical Summary Report #95-04, October 1994,
    <http://ftp.cs.wisc.edu/Approx/symmetries.pdf>.
    """

    def __init__(self, r, alpha):
        assert r ** 2 >= 1 / 3

        self.name = "Waldron"
        self.degree = 3

        R = r / numpy.sqrt(3 * r ** 2 - 1)

        beta = alpha + numpy.pi / 2

        self.points = numpy.array(
            [
                [+r * numpy.cos(alpha), +r * numpy.sin(alpha)],
                [-r * numpy.cos(alpha), -r * numpy.sin(alpha)],
                [+R * numpy.cos(beta), +R * numpy.sin(beta)],
                [-R * numpy.cos(beta), -R * numpy.sin(beta)],
            ]
        )

        self.weights = numpy.array(
            [2 / 3 / r ** 2, 2 / 3 / r ** 2, 2 / 3 / R ** 2, 2 / 3 / R ** 2]
        )
        return
