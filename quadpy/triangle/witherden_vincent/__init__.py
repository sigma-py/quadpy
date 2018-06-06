# -*- coding: utf-8 -*-
#
import json
import os

from ..helpers import untangle2


class WitherdenVincent(object):
    """
    F.D. Witherden, P.E. Vincent,
    On the identification of symmetric quadrature rules for finite
    element methods,
    Computers & Mathematics with Applications,
    Volume 69, Issue 10, May 2015, Pages 1232–1241,
    <https://doi.org/10.1016/j.camwa.2015.03.017>.

    Abstract:
    In this paper we describe a methodology for the identification of symmetric
    quadrature rules inside of quadrilaterals, triangles, tetrahedra, prisms,
    pyramids, and hexahedra. The methodology is free from manual intervention
    and is capable of identifying a set of rules with a given strength and a
    given number of points. We also present polyquad which is an implementation
    of our methodology. Using polyquad v1.0 we proceed to derive a complete set
    of symmetric rules on the aforementioned domains. All rules possess purely
    positive weights and have all points inside the domain. Many of the rules
    appear to be new, and an improvement over those tabulated in the
    literature.
    """

    def __init__(self, degree, symbolic=False):
        self.name = "WV({})".format(degree)

        this_dir = os.path.dirname(os.path.realpath(__file__))
        filename = "wv{:02d}.json".format(degree)
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)

        self.degree = data.pop("degree")

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
