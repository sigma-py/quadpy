import json

import numpy

from ..helpers import expand_symmetries, plot_disks
from ..tn import TnScheme, get_vol, transform

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class T2Scheme(TnScheme):
    def __init__(
        self,
        name,
        symmetry_data,
        degree,
        source=None,
        tol=1.0e-14,
        comments=None,
        weight_factor=None,
    ):
        self.symmetry_data = symmetry_data

        points, weights = expand_symmetries(symmetry_data)
        if weight_factor is not None:
            weights *= weight_factor

        assert points.shape[0] == 3, f"{name}, {points.shape}"
        super().__init__(name, 2, weights, points, degree, source, tol, comments)
        self.domain = "T2"

    def plot(
        self,
        triangle=numpy.array([[-0.5, 0.0], [+0.5, 0.0], [0, 0.5 * (numpy.sqrt(3))]]),
        show_axes=False,
    ):
        """Shows the quadrature points on a given triangle. The size of the circles
        around the points coincides with their weights.
        """
        from matplotlib import pyplot as plt

        plt.plot(triangle[:, 0], triangle[:, 1], "-k")
        plt.plot(
            [triangle[-1, 0], triangle[0, 0]], [triangle[-1, 1], triangle[0, 1]], "-k"
        )

        if not show_axes:
            plt.gca().set_axis_off()

        transformed_pts = transform(self.points, triangle.T).T

        vol = get_vol(triangle)
        plot_disks(plt, transformed_pts, self.weights, vol)

        plt.axis("equal")


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]
    d = content["data"]
    weight_factor = content["weight factor"] if "weight factor" in content else None
    return T2Scheme(name, d, degree, source, tol, weight_factor=weight_factor)


def _scheme_from_dict(content, source=None):
    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return T2Scheme(
        content["name"],
        weights,
        points,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )
