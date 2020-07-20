import json
import warnings

import numpy
import sympy

from ..helpers import plot_disks
from ..tn import TnScheme, get_vol, transform


class T2Scheme(TnScheme):
    def __init__(self, name, weights, points, degree, source=None, tol=1.0e-14):
        super().__init__(name, 2, weights, points, degree, source, tol)
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

        transformed_pts = transform(self.points.T, triangle.T).T

        vol = get_vol(triangle)
        plot_disks(plt, transformed_pts, self.weights, vol)

        plt.axis("equal")


def _s3(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return numpy.full((1, 3), frac(1, 3))


def _s21(a):
    a = numpy.array(a)
    b = 1 - 2 * a
    return numpy.array([[a, a, b], [a, b, a], [b, a, a]])


def _s111ab(a, b):
    c = 1 - a - b
    out = numpy.array(
        [[a, b, c], [c, a, b], [b, c, a], [b, a, c], [c, b, a], [a, c, b]]
    )
    out = numpy.swapaxes(out, 0, 1)
    return out


def _rot_ab(a, b):
    c = 1 - a - b
    out = numpy.array([[a, b, c], [c, a, b], [b, c, a]])
    out = numpy.swapaxes(out, 0, 1)
    return out


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)


def untangle2(data, symbolic=False):
    points = []
    weights = []

    if "s3" in data:
        d = numpy.array(data["s3"]).T
        points.append(_s3(symbolic).T)
        weights.append(numpy.tile(d[0], 1))

    if "s2" in data:
        d = numpy.array(data["s2"]).T
        s2_data = _s21(d[1])
        points.append(_collapse0(s2_data))
        weights.append(numpy.tile(d[0], 3))

    if "s1" in data:
        d = numpy.array(data["s1"]).T
        s1_data = _s111ab(*d[1:])
        points.append(_collapse0(s1_data))
        weights.append(numpy.tile(d[0], 6))

    if "rot" in data:
        d = numpy.array(data["rot"]).T
        rot_data = _rot_ab(*d[1:])
        points.append(_collapse0(rot_data))
        weights.append(numpy.tile(d[0], 3))

    points = numpy.column_stack(points).T
    weights = numpy.concatenate(weights)
    return points, weights


def s3(weight):
    symbolic = isinstance(weight, float)
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return numpy.array([weight]), numpy.full((1, 3), frac(1, 3))


def s2(*data):
    w, a = numpy.array(data).T
    b = 1 - 2 * a
    points = _stack_first_last([[a, a, b], [a, b, a], [b, a, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def s1(*data):
    w, a, b = numpy.array(data).T
    c = 1 - a - b
    points = _stack_first_last(
        [[a, b, c], [c, a, b], [b, c, a], [b, a, c], [c, b, a], [a, c, b]]
    )
    weights = numpy.tile(w, 6)
    return weights, points


def r(*data):
    w, r = numpy.array(data).T
    a = r + (1 - r) / 3
    b = (1 - a) / 2
    points = _stack_first_last([[a, b, b], [b, a, b], [b, b, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def rot_ab(*data):
    w, a, b = numpy.array(data).T
    c = 1 - a - b
    points = _stack_first_last([[a, b, c], [c, a, b], [b, c, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def mirror(*data):
    w, a, b = numpy.array(data).T
    c = 1 - a - b
    points = _stack_first_last([[a, b, c], [b, a, c]])
    weights = numpy.tile(w, 2)
    return weights, points


def alpha(*data):
    """From the article Liu-Vinokur:

    mu_i = (1 + (n-1) alpha) / n,
    mu_j = (1 - alpha) / n    for j!=i,

    where n is the number of vertices
    """
    w, alpha = numpy.array(data).T
    a = (1 + 2 * alpha) / 3
    b = (1 - alpha) / 3
    points = _stack_first_last([[a, b, b], [b, a, b], [b, b, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def gamma_delta(*data):
    """From the article Liu-Vinokur:

    mu_i = (1 + (n-1) gamma - delta) / n,
    mu_j = (1 + (n-1) delta - gamma) / n,
    mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

    where n is the number of vertices
    """
    w, gamma, delta = numpy.array(data).T
    a = (1 + 2 * gamma - delta) / 3
    b = (1 + 2 * delta - gamma) / 3
    c = (1 - gamma - delta) / 3
    points = _stack_first_last(
        [[a, b, c], [c, a, b], [b, c, a], [a, c, b], [b, a, c], [c, b, a]]
    )
    weights = numpy.tile(w, 6)
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


def _s3_alt(data):
    return numpy.full((3, 1), 1 / 3)


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


def _rot_ab_alt(data):
    a, b = data
    c = 1 - a - b
    points = numpy.array([[a, b, c], [c, a, b], [b, c, a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {"s1": _s1, "s2": _s2, "s3": _s3_alt, "rot_ab": _rot_ab_alt}[key]
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
    if "test_tolerance" in content:
        tol = content["test_tolerance"]
        if tol > 1.0e-12:
            warnings.warn(f"The {name} scheme has low precision ({tol:.3e}).")
    else:
        tol = 1.0e-14

    points, weights = expand_symmetries(content["data"])

    if weight_factor is not None:
        weights *= weight_factor

    return T2Scheme(name, weights, points, degree, source, tol)
