# -*- coding: utf-8 -*-
#
import json
import os

from ..helpers import untangle2, cartesian_to_spherical


class Lebedev(object):
    """
    Sphere integration schemes from

    Lebedev, V. I. (1976),
    Quadratures on a sphere,
    Zh. Vȳchisl. Mat. Mat. Fiz. 16 (2): 293–306,
    <https://doi.org/10.1016/0041-5553(76)90100-2>.

    <https://en.wikipedia.org/wiki/Lebedev_quadrature>
    <https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html>
    """

    # It's a little unclear how to best store the original data. By Burkhardt,
    # it is given in terms of the spherical coordinates phi and theta, however
    # those angles are not well suited to express the symmetry. For that, this
    # code converts the spherical coordinates into Cartesians, applies the
    # symmetry transformations, and converts back.
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
