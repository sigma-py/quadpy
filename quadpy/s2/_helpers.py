import json
import warnings

import numpy

from ..helpers import QuadratureScheme, expand_symmetries, plot_disks


class S2Scheme(QuadratureScheme):
    def __init__(
        self, name, symmetry_data, degree: int, source=None, tol=1.0e-14, weight_factor=None
    ):
        self.domain = "S2"
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(d)
        assert points.shape[0] == 2
        assert points.shape[1] == weights.shape[0], f"{points.shape}, {weights.shape}"
        if weight_factor is not None:
            weights *= weight_factor
        super().__init__(name, weights, points, degree, source, tol)

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        ax = plt.gca()
        # change default range so that new disks will work
        plt.axis("equal")
        ax.set_xlim((-1.5, 1.5))
        ax.set_ylim((-1.5, 1.5))

        if not show_axes:
            ax.set_axis_off()

        disk1 = plt.Circle((0, 0), 1, color="k", fill=False)
        ax.add_artist(disk1)

        plot_disks(plt, self.points.T, self.weights, numpy.pi)
        return

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points.T)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.array(f((rr + center).T))
        return numpy.pi * numpy.array(radius) ** 2 * dot(ff, self.weights)


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    if "_ERR" in content:
        warnings.warn(content["_ERR"])

    weight_factor = content["weight factor"] if "weight factor" in content else None

    return S2Scheme(
        name, content["data"], degree, source, tol, weight_factor=weight_factor
    )


def _scheme_from_dict(content, source=None):
    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return S2Scheme(
        content["name"],
        weights,
        points,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )
