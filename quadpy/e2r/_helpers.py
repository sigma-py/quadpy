import json
import math

import numpy

from ..helpers import QuadratureScheme, expand_symmetries, plot_disks


class E2rScheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        self.domain = "E2r"
        assert points.shape[0] == 2
        super().__init__(name, weights, points, degree, source, tol)

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        ax = plt.gca()
        plt.axis("equal")

        if not show_axes:
            ax.set_axis_off()

        I0 = 2 * math.pi

        plot_disks(plt, self.points.T, self.weights, I0)

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = 2 * math.pi
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return E2rScheme(name, weights, points, degree, source, tol)
