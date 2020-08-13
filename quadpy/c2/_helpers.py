import json

import numpy

from ..cn import CnScheme
from ..cn import ncube_points as rectangle_points
from ..cn import transform
from ..helpers import expand_symmetries, plot_disks
from ..tn import get_vol


class C2Scheme(CnScheme):
    def __init__(
        self, name, weights, points, degree, source=None, tol=1.0e-14, comments=None
    ):
        super().__init__(name, 2, weights, points, degree, source, tol, comments)
        self.domain = "C2"

    def plot(self, quad=rectangle_points([0.0, 1.0], [0.0, 1.0]), show_axes=False):
        """Shows the quadrature points on a given quad. The area of the disks
        around the points coincides with their weights.
        """
        from matplotlib import pyplot as plt

        def plot_segment(a, b):
            plt.plot((a[0], b[0]), (a[1], b[1]), "-k")

        plot_segment(quad[0][0], quad[1][0])
        plot_segment(quad[1][0], quad[1][1])
        plot_segment(quad[1][1], quad[0][1])
        plot_segment(quad[0][1], quad[0][0])

        if not show_axes:
            plt.gca().set_axis_off()

        transformed_pts = transform(self.points, quad)

        # compute volume by splitting it in two triangles
        vol = get_vol(numpy.array([quad[0][0], quad[1][0], quad[0][1]])) + get_vol(
            numpy.array([quad[0][0], quad[0][1], quad[1][1]])
        )
        plot_disks(plt, transformed_pts, self.weights, vol)
        plt.axis("equal")
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return C2Scheme(name, weights, points, degree, source, tol)


def _scheme_from_dict(content, source=None):
    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return C2Scheme(
        content["name"],
        weights,
        points,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )
