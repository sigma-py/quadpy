import json

import numpy

from ..helpers import plot_disks
from ..tn import TnScheme, get_vol, transform


class T2Scheme(TnScheme):
    def __init__(
        self, name, weights, points, degree, source=None, tol=1.0e-14, comments=None
    ):
        self.domain = "T2"
        assert points.shape[0] == 3, f"{name}, {points.shape}"
        super().__init__(name, 2, weights, points, degree, source, tol, comments)

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


def _s3(data):
    return numpy.full((3, 1), 1 / 3)


def _vertex(data):
    return numpy.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


def _s2(a):
    a = numpy.array(a)
    b = 1 - 2 * a
    return numpy.array([[a, a, b], [a, b, a], [b, a, a]])


def _s1(data):
    a, b = numpy.asarray(data)
    c = 1 - a - b
    points = numpy.array(
        [[a, b, c], [c, a, b], [b, c, a], [b, a, c], [c, b, a], [a, c, b]]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _rot_ab(data):
    a, b = data
    c = 1 - a - b
    points = numpy.array([[a, b, c], [c, a, b], [b, c, a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _swap_ab(data):
    a, b = data
    c = 1 - a - b
    points = numpy.array([[a, b, c], [b, a, c]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s2_static(a):
    a = numpy.asarray(a)
    b = 1 - 2 * a
    points = numpy.array([[a, a, b]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "s1": _s1,
            "s2": _s2,
            "s3": _s3,
            "rot": _rot_ab,
            "rot_ab": _rot_ab,
            "swap_ab": _swap_ab,
            "s2_static": _s2_static,
            "vertex": _vertex,
        }[key]
        pts = fun(numpy.asarray(points_raw))

        counts.append(pts.shape[1])
        pts = pts.reshape(pts.shape[0], -1)
        points.append(pts)

    points = numpy.ascontiguousarray(numpy.concatenate(points, axis=1))
    return points, counts


def expand_symmetries(data):
    # separate points and weights
    points_raw = {}
    weights_raw = []
    for key, values in data.items():
        weights_raw.append(values[0])
        points_raw[key] = values[1:]

    points, counts = expand_symmetries_points_only(points_raw)
    weights = numpy.concatenate(
        [numpy.tile(values, count) for count, values in zip(counts, weights_raw)]
    )
    return points, weights


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return T2Scheme(name, weights, points, degree, source, tol)


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
