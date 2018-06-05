# -*- coding: utf-8 -*-
#
import json
import os

from ..helpers import untangle2


class Papanicolopulos(object):
    """
    Stefanos-Aldo Papanicolopulos,
    New fully symmetric and rotationally symmetric cubature rules on the
    triangle using minimal orthonormal bases,
    Journal of Computational and Applied Mathematics,
    Volume 294, 1 March 2016, Pages 39–48,
    <https://doi.org/10.1016/j.cam.2015.08.001>,
    <https://arxiv.org/abs/1411.5631>.
    """

    def __init__(self, variant, index):
        self.name = "Papanicolopulos({}, {})".format(variant, index)

        this_dir = os.path.dirname(os.path.realpath(__file__))
        if variant == "fs":
            filename = "full{:02d}.json".format(index)
        else:
            assert variant == "rot"
            # ERR the first 8 schemes are flawed by round-off error
            assert index >= 8
            filename = "rot{:02d}.json".format(index)

        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)
        self.degree = data.pop("degree")

        self.data = data
        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
