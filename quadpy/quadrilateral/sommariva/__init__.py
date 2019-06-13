# -*- coding: utf-8 -*-
#
"""
Quadrilateral quadrature by Alvise Sommariva,
<http://www.math.unipd.it/~alvise/sets.html>.
"""
import json
import os

import numpy

from ..helpers import QuadrilateralScheme


def sommariva(index):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = "sommariva_{:02d}.json".format(index)
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    data = numpy.array(data["data"])
    points = data[:, :2]
    weights = data[:, 2]
    return QuadrilateralScheme("Sommariva {}".format(index), degree, weights, points)


Sommariva = {k: lambda: sommariva(k) for k in range(1, 56)}
