import json
import math

import numpy

from ..helpers import QuadratureScheme, plot_disks


class E2r2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        self.domain = "E2r2"
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
        ref_vol = math.pi
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))


def _zero(data):
    return numpy.array([[0.0], [0.0]])


def _pm2(data):
    x, y = data
    points = numpy.array([[+x, +y], [+x, -y], [-x, +y], [-x, -y]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmx(r):
    zero = numpy.zeros_like(r)
    points = numpy.array([[+r, zero], [-r, zero]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmy(r):
    zero = numpy.zeros_like(r)
    points = numpy.array([[zero, +r], [zero, -r]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s8(data):
    a, b = data
    points = numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s4(a):
    points = numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s40(a):
    zero = numpy.zeros_like(a)
    points = numpy.array([[+a, zero], [-a, zero], [zero, +a], [zero, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "zero": _zero,
            "pm2": _pm2,
            "pmx": _pmx,
            "pmy": _pmy,
            "s40": _s40,
            "s4": _s4,
            "s8": _s8,
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

    return E2r2Scheme(name, weights, points, degree, source, tol)
