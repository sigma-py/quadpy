# -*- coding: utf-8 -*-
#
import json
import os

from ..helpers import untangle2


class WandzuraXiao(object):
    """
    S. Wandzura, H. Xiao,
    Symmetric quadrature rules on a triangle,
    Computers & Mathematics with Applications,
    Volume 45, Issue 12, June 2003, Pages 1829-1840,
    <https://doi.org/10.1016/S0898-1221(03)90004-6>.

    Abstract:
    We present a class of quadrature rules on triangles in R2 which, somewhat
    similar to Gaussian rules on intervals in R1, have rapid convergence,
    positive weights, and symmetry. By a scheme combining simple group theory
    and numerical optimization, we obtain quadrature rules of this kind up to
    the order 30 on triangles. This scheme, essentially a formalization and
    generalization of the approach used by Lyness and Jespersen over 25 years
    ago, can be easily extended to other regions in R2 and surfaces in higher
    dimensions, such as squares, spheres. We present example formulae and
    relevant numerical results.

    Note that in the above article, the authors present the coordinates in the
    symmetric triangle [[-0.5, -sqrt(3)/2], [-0.5, +sqrt(3)/2], [1, 0]]. These
    have been transformed to barycentric coordinates here.
    """

    def __init__(self, index):
        self.name = "WX({})".format(index)

        this_dir = os.path.dirname(os.path.realpath(__file__))
        filename = "wx{:02d}.json".format(index)
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)
        self.degree = data.pop("degree")

        self.data = data
        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
