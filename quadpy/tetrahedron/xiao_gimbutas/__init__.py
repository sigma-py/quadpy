# -*- coding: utf-8 -*-
#
import json
import os

from ..helpers import untangle2


class XiaoGimbutas(object):
    """
    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature rules in
    two and higher dimensions,
    Computers & Mathematics with Applications,
    Volume 59, Issue 2, January 2010, Pages 663–676,
    <https://doi.org/10.1016/j.camwa.2009.10.027>.

    Abstract:
    We present a numerical algorithm for the construction of efficient,
    high-order quadratures in two and higher dimensions. Quadrature rules
    constructed via this algorithm possess positive weights and interior nodes,
    resembling the Gaussian quadratures in one dimension. In addition, rules
    can be generated with varying degrees of symmetry, adaptable to individual
    domains. We illustrate the performance of our method with numerical
    examples, and report quadrature rules for polynomials on triangles,
    squares, and cubes, up to degree 50. These formulae are near optimal in the
    number of nodes used, and many of them appear to be new.

    Data adapted from
    <https://people.sc.fsu.edu/~jburkardt/f_src/triangle_symq_rule/triangle_symq_rule.f90>.
    """

    def __init__(self, degree, symbolic=False):
        self.name = "XG({})".format(degree)

        this_dir = os.path.dirname(os.path.realpath(__file__))
        filename = "xg{:02d}.json".format(degree)
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)

        self.degree = data.pop("degree")

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
