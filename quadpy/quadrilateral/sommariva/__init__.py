# -*- coding: utf-8 -*-
#
import json
import os
import re

import numpy


class Sommariva(object):
    """
    """

    def __init__(self, index):
        self.name = "Sommariva({})".format(index)

        this_dir = os.path.dirname(os.path.realpath(__file__))

        m = re.match("([0-9]+)([a-z]*)", index)
        filename = "sommariva_{:02d}{}.json".format(int(m.group(1)), m.group(2))
        with open(os.path.join(this_dir, filename), "r") as f:
            data = json.load(f)

        self.degree = data.pop("degree")

        data = numpy.array(data["data"])
        self.points = data[:, :2]
        self.weights = data[:, 2]
        return
