import json
import warnings

import numpy

from ..helpers import QuadratureScheme, plot_disks


class S2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree: int, source=None, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "S2"

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


def _z():
    return numpy.array([[0.0, 0.0]])


def _s8(a, b):
    return numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )


def _s4(a):
    return numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])


def _s40(a):
    return numpy.array([[+a, 0.0], [-a, 0.0], [0.0, +a], [0.0, -a]])


def _pm(a, b):
    return numpy.array([[+a, +b], [-a, +b], [+a, -b], [-a, -b]])


def _pmx(x):
    return numpy.array([[+x, 0], [-x, 0]])


def _pmy(y):
    return numpy.array([[0, +y], [0, -y]])


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


def _zero_alt(data):
    return numpy.array([[0.0], [0.0]])


def _fsd_alt(data):
    a, b = numpy.asarray(data)
    points = numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "pm": _pm_alt,
            "pmx": _pmx_alt,
            "zero": _zero_alt,
            "fsd": _fsd_alt
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

    # TODO remove this once points are expected as points.T in all functions
    points = points.T
    return points, weights


def _read(filepath, source, weight_factor=None):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    if tol > 1.0e-12:
        warnings.warn(f"The {name} scheme has low precision ({tol:.3e}).")

    points, weights = expand_symmetries(content["data"])

    if weight_factor is not None:
        weights *= weight_factor

    return S2Scheme(name, weights, points, degree, source, tol)
