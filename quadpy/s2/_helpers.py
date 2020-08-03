import json
import warnings

import numpy
import sympy

from ..helpers import QuadratureScheme, plot_disks


class S2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree: int, source=None, tol=1.0e-14):
        self.domain = "S2"
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

        plot_disks(plt, self.points, self.weights, numpy.pi)
        return

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.array(f((rr + center).T))
        return numpy.pi * numpy.array(radius) ** 2 * dot(ff, self.weights)


def _pma_alt(data):
    a = numpy.asarray(data)
    points = numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pm_alt(data):
    a, b = numpy.asarray(data)
    points = numpy.array([[+a, +b], [-a, +b], [+a, -b], [-a, -b]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmx_alt(data):
    a = numpy.asarray(data)
    zero = numpy.zeros_like(a)
    points = numpy.array([[+a, zero], [-a, zero]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmy_alt(data):
    a = numpy.asarray(data)
    zero = numpy.zeros_like(a)
    points = numpy.array([[zero, +a], [zero, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _zero_alt(data):
    if data.dtype == sympy.Basic:
        return numpy.array([[0], [0]])
    return numpy.array([[0.0], [0.0]])


def _fsd_alt(data):
    a, b = numpy.asarray(data)
    points = numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s40_alt(data):
    a = numpy.asarray(data)
    zero = numpy.zeros_like(a)
    points = numpy.array([[+a, zero], [-a, zero], [zero, +a], [zero, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "pm": _pm_alt,
            "pmx": _pmx_alt,
            "pmy": _pmy_alt,
            "zero": _zero_alt,
            "fsd": _fsd_alt,
            "s40": _s40_alt,
            "pma": _pma_alt,
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
        values = numpy.asarray(values)
        weights_raw.append(values[0])
        points_raw[key] = values[1:]

    points, counts = expand_symmetries_points_only(points_raw)
    weights = numpy.concatenate(
        [numpy.tile(values, count) for count, values in zip(counts, weights_raw)]
    )

    # TODO remove this once points are expected as points.T in all functions
    points = points.T
    return points, weights


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    if "_ERR" in content:
        warnings.warn(content["_ERR"])

    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return S2Scheme(name, weights, points, degree, source, tol)
