import json
import math

import numpy

from ..helpers import QuadratureScheme, expand_symmetries, plot_disks


class E2r2Scheme(QuadratureScheme):
    def __init__(
        self, name, symmetry_data, degree, source, tol=1.0e-14, weight_factor=None
    ):
        self.domain = "E2r2"
        points, weights = expand_symmetries(symmetry_data)
        assert points.shape[0] == 2
        if weight_factor is not None:
            weights *= weight_factor
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
        ref_vol = math.pi
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]
    data = content["data"]

    weight_factor = content["weight factor"] if "weight factor" in content else None

    return E2r2Scheme(
        name, data, degree, source, tol, weight_factor=weight_factor
    )
