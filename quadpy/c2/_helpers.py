import json

import numpy
import orthopy

from .. import helpers
from ..cn import CnScheme
from ..cn import ncube_points as rectangle_points
from ..cn import transform
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
        helpers.plot_disks(plt, transformed_pts, self.weights, vol)
        plt.axis("equal")
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)

    def compute_residuals(self, level):
        evaluator = orthopy.cn.Eval(self.points)

        quad = rectangle_points([-1.0, +1.0], [-1.0, +1.0])

        max_res = []
        for k in range(level + 1):
            vals = next(evaluator)
            approximate = self.integrate(lambda x: vals, quad)
            exact = 2.0 if k == 0 else 0.0
            res = numpy.abs(approximate - exact)
            max_res += [numpy.max(res)]

        return numpy.array(max_res)


def zero(weight):
    return numpy.array([weight]), numpy.array([[0, 0]])


def pmx(*data):
    w, x, y = numpy.array(data).T
    points = _stack_first_last([[+x, y], [-x, y]])
    weights = numpy.tile(w, 2)
    return weights, points


def pm(*data):
    w, s, t = numpy.array(data).T
    points = _stack_first_last([[+s, +t], [-s, -t]])
    weights = numpy.tile(w, 2)
    return weights, points


def pm2(*data):
    w, x, y = numpy.array(data).T
    points = _stack_first_last([[+x, +y], [+x, -y], [-x, +y], [-x, -y]])
    weights = numpy.tile(w, 4)
    return weights, points


def symm_r0(*data):
    data = numpy.array(data)
    w, r = data.T
    zero = numpy.zeros(w.shape[0], dtype=r.dtype)
    points = _stack_first_last([[+r, zero], [-r, zero], [zero, +r], [zero, -r]])
    weights = numpy.tile(w, 4)
    return weights, points


def symm_s(*data):
    w, s = numpy.array(data).T
    points = _stack_first_last([[+s, +s], [+s, -s], [-s, +s], [-s, -s]])
    weights = numpy.tile(w, 4)
    return weights, points


def symm_s_t(*data):
    w, s, t = numpy.array(data).T
    points = _stack_first_last(
        [[+s, +t], [-s, +t], [+s, -t], [-s, -t], [+t, +s], [-t, +s], [+t, -s], [-t, -s]]
    )
    weights = numpy.tile(w, 8)
    return weights, points


def _stack_first_last(arr):
    """Stacks an input array of shape (i, j, k) such that the output array is of shape
    (i*k, j).
    """
    arr = numpy.swapaxes(arr, 0, 1)
    return arr.reshape(arr.shape[0], -1).T


def concat(*data):
    weights = numpy.concatenate([t[0] for t in data])
    points = numpy.vstack([t[1] for t in data])
    return weights, points


def _symm_s_t(data):
    s, t = numpy.array(data)
    points = numpy.array(
        [[+s, +t], [-s, +t], [+s, -t], [-s, -t], [+t, +s], [-t, +s], [+t, -s], [-t, -s]]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_s(a):
    points = numpy.array([[+a, +a], [+a, -a], [-a, +a], [-a, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_r0(r):
    zero = numpy.zeros_like(r)
    points = numpy.array([[+r, zero], [-r, zero], [zero, +r], [zero, -r]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s4(data):
    a, b = data
    points = numpy.array([[+a, +b], [-a, -b], [-b, +a], [+b, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _zero(data):
    return numpy.array([[0.0], [0.0]])


def _pm2(data):
    x, y = data
    points = numpy.array([[+x, +y], [+x, -y], [-x, +y], [-x, -y]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pm(data):
    x, y = data
    points = numpy.array([[+x, +y], [-x, -y]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmx2(data):
    x, y = data
    points = numpy.array([[+x, y], [-x, y]])
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


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "symm_s_t": _symm_s_t,
            "symm_s": _symm_s,
            "symm_r0": _symm_r0,
            "s4": _s4,
            "zero": _zero,
            "pm2": _pm2,
            "pm": _pm,
            "pmx": _pmx,
            "pmx2": _pmx2,
            "pmy": _pmy,
            "plain": lambda vals: vals.reshape(2, 1, -1),
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
