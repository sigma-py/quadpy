# -*- coding: utf-8 -*-
#
import json
import os

from ..helpers import untangle2, cartesian_to_spherical


class Lebedev(object):
    """
    Sphere integration schemes from a series of publications, in chronological order:

    V.I. Lebedev,
    Values of the nodes and weights of ninth to seventeenth order Gauss-Markov
    quadrature formulae invariant under the octahedron group with inversion,
    Computational Mathematics and Mathematical Physics, Vol. 15, 1975, pp. 44-51.

    Lebedev, V.I. (1976),
    Quadratures on a sphere,
    Zh. Vȳchisl. Mat. Mat. Fiz. 16 (2): 293–306,
    <https://doi.org/10.1016/0041-5553(76)90100-2>.

    V.I. Lebedev,
    Spherical quadrature formulas exact to orders 25-29,
    Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107.

    V.I. Lebedev and A.L. Skorokhodov,
    Quadrature formulas of orders 41, 47, and 53 for the sphere,
    Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592.

    V.I. Lebedev,
    A quadrature formula for the sphere of 59th algebraic order of accuracy,
    Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286.

    V.I. Lebedev and D.N. Laikov,
    A quadrature formula for the sphere of the 131st algebraic order of accuracy,
    Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.

    <https://en.wikipedia.org/wiki/Lebedev_quadrature>
    <https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html>
    """
    def __init__(self, degree):
        self.name = "Lebedev({})".format(degree)

        this_dir = os.path.dirname(os.path.realpath(__file__))
        filename = "lebedev_{:03d}.json".format(degree)
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)

        self.degree = data.pop("degree")

        self.points, self.weights = untangle2(data)
        self.azimuthal_polar = cartesian_to_spherical(self.points)
        return
