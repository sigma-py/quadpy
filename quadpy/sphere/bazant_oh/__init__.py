# -*- coding: utf-8 -*-
#
import json
import os
import re

from ..helpers import cartesian_to_spherical, untangle2


class BazantOh(object):
    """
    P. Bažant, B.H. Oh,
    Efficient Numerical Integration on the Surface of a Sphere,
    ZAMM, Volume 66, Issue 1, 1986, Pages 37-49,
    https://doi.org/10.1002/zamm.19860660108.

    Abstract:
    Several new numerical integration formulas on the surface of a sphere in three
    dimensions are derived. These formulas are superior to the existing ones in that for
    the same degree of approximation they require fewer integration points for functions
    with central or planar symmetry. Furthermore, a general method of deriving the
    integration formulas, which achieves conceptual simplicity at the expense of
    extensive numerical work left for a computer, is demonstrated. In this method, the
    coefficients of the integration formula are determined from a system of linear
    algebraic equations directly representing the conditions for a certain number of
    terms of the three‐dimensional Taylor series expansion of the integrated function
    about the center of the sphere to vanish, while the unknown locations of the
    integration points are determined from a condition for the next term (or terms) of
    the expansion to vanish, and if it cannot be made to vanish, then from a condition
    for minimizing the magnitude of this term (or these terms). Finally, we formulate a
    new condition of optimality of the integration formulas which is important for the
    integration error in certain applications.
    """

    def __init__(self, index):
        self.name = "BazantOh({})".format(index)

        this_dir = os.path.dirname(os.path.realpath(__file__))

        m = re.match("([0-9]+)([a-z]*)", index)
        filename = "bazant_oh_{:03d}{}.json".format(int(m.group(1)), m.group(2))
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)

        self.degree = data.pop("degree")

        self.points, self.weights = untangle2(data)
        self.azimuthal_polar = cartesian_to_spherical(self.points)
        return
