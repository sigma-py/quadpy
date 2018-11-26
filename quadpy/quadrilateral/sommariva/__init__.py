# -*- coding: utf-8 -*-
#
import json
import os

import numpy


class Sommariva(object):
    """
    Quad quadrature by Alvise Sommariva, <http://www.math.unipd.it/~alvise/sets.html>.
    """

    def __init__(self, index):
        self.name = "Sommariva({})".format(index)

        this_dir = os.path.dirname(os.path.realpath(__file__))

        filename = "sommariva_{:02d}.json".format(index)
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)

        self.degree = data.pop("degree")

        data = numpy.array(data["data"])
        self.points = data[:, :2]
        self.weights = data[:, 2]
        return
