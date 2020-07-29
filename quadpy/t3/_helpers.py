import json

import numpy
import sympy

from ..helpers import backend_to_function
from ..tn import TnScheme, get_vol, transform


class T3Scheme(TnScheme):
    def __init__(self, name, weights, points, degree, source=None, tol=1.0e-14):
        self.domain = "T3"
        super().__init__(name, 2, weights, points, degree, source, tol)

    def show(
        self,
        tet=numpy.array(
            [
                [+1, 0, -1.0 / numpy.sqrt(2.0)],
                [-1, 0, -1.0 / numpy.sqrt(2.0)],
                [0, +1, +1.0 / numpy.sqrt(2.0)],
                [0, -1, +1.0 / numpy.sqrt(2.0)],
            ]
        ),
        backend="vtk",
        render=True,
    ):
        edges = numpy.array([[tet[i], tet[j]] for i in range(4) for j in range(i)])
        edges = numpy.moveaxis(edges, 1, 2)
        backend_to_function[backend](
            transform(self.points.T, tet.T).T,
            self.weights,
            get_vol(tet),
            edges,
            render=render,
        )
        return


def _s4(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return numpy.full((1, 4), frac(1, 4))


def _s31(a):
    b = 1 - 3 * a
    return numpy.array([[a, a, a, b], [a, a, b, a], [a, b, a, a], [b, a, a, a]])


def _s22(a):
    b = 0.5 - a
    return numpy.array(
        [
            [a, a, b, b],
            [a, b, a, b],
            [b, a, a, b],
            [a, b, b, a],
            [b, a, b, a],
            [b, b, a, a],
        ]
    )


def _s211(a, b):
    c = 1.0 - 2 * a - b
    return numpy.array(
        [
            [a, a, b, c],
            [a, b, a, c],
            [b, a, a, c],
            [a, b, c, a],
            [b, a, c, a],
            [b, c, a, a],
            [a, a, c, b],
            [a, c, a, b],
            [c, a, a, b],
            [a, c, b, a],
            [c, a, b, a],
            [c, b, a, a],
        ]
    )


def _s1111(a, b, c):
    d = 1.0 - a - b - c
    return numpy.array(
        [
            [a, b, c, d],
            [a, b, d, c],
            [a, c, b, d],
            [a, c, d, b],
            [a, d, b, c],
            [a, d, c, b],
            [b, a, c, d],
            [b, a, d, c],
            [b, c, a, d],
            [b, c, d, a],
            [b, d, a, c],
            [b, d, c, a],
            [c, a, b, d],
            [c, a, d, b],
            [c, b, a, d],
            [c, b, d, a],
            [c, d, a, b],
            [c, d, b, a],
            [d, a, b, c],
            [d, a, c, b],
            [d, b, a, c],
            [d, b, c, a],
            [d, c, a, b],
            [d, c, b, a],
        ]
    )


def untangle2(data):
    points = []
    weights = []

    if "s4" in data:
        assert len(data["s4"]) == 1
        w = numpy.array(data["s4"]).T
        points.append(_s4(symbolic=False))
        weights.append(w[0])

    if "s31" in data:
        d = numpy.array(data["s31"]).T
        s31_data = numpy.moveaxis(_s31(d[1]), 0, 1)
        points.append(_collapse0(s31_data).T)
        weights.append(numpy.tile(d[0], 4))

    if "s22" in data:
        d = numpy.array(data["s22"]).T
        s22_data = numpy.moveaxis(_s22(d[1]), 0, 1)
        points.append(_collapse0(s22_data).T)
        weights.append(numpy.tile(d[0], 6))

    if "s211" in data:
        d = numpy.array(data["s211"]).T
        s211_data = numpy.moveaxis(_s211(*d[1:]), 0, 1)
        points.append(_collapse0(s211_data).T)
        weights.append(numpy.tile(d[0], 12))

    if "s1111" in data:
        d = numpy.array(data["s1111"]).T
        s1111_data = numpy.moveaxis(_s1111(*d[1:]), 0, 1)
        points.append(_collapse0(s1111_data).T)
        weights.append(numpy.tile(d[0], 24))

    points = numpy.concatenate(points)
    weights = numpy.concatenate(weights)

    return points, weights


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)


def s4(weight):
    symbolic = not isinstance(weight, float)
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return numpy.array([weight]), numpy.full((1, 4), frac(1, 4))


def s31(*data):
    w, a = numpy.array(data).T
    b = 1 - 3 * a
    points = _stack_first_last([[a, a, a, b], [a, a, b, a], [a, b, a, a], [b, a, a, a]])
    weights = numpy.tile(w, 4)
    return weights, points


def s22(*data):
    w, a = numpy.array(data).T
    b = (1 - 2 * a) / 2
    points = _stack_first_last(
        [
            [a, a, b, b],
            [a, b, a, b],
            [b, a, a, b],
            [a, b, b, a],
            [b, a, b, a],
            [b, b, a, a],
        ]
    )
    weights = numpy.tile(w, 6)
    return weights, points


def s211(*data):
    w, a, b = numpy.array(data).T
    c = 1 - 2 * a - b
    points = _stack_first_last(
        [
            [a, a, b, c],
            [a, b, a, c],
            [b, a, a, c],
            [a, b, c, a],
            [b, a, c, a],
            [b, c, a, a],
            [a, a, c, b],
            [a, c, a, b],
            [c, a, a, b],
            [a, c, b, a],
            [c, a, b, a],
            [c, b, a, a],
        ]
    )
    weights = numpy.tile(w, 12)
    return weights, points


def s1111(*data):
    w, a, b, c = numpy.array(data).T
    d = 1 - a - b - c
    points = _stack_first_last(
        [
            [a, b, c, d],
            [a, b, d, c],
            [a, c, b, d],
            [a, c, d, b],
            [a, d, b, c],
            [a, d, c, b],
            [b, a, c, d],
            [b, a, d, c],
            [b, c, a, d],
            [b, c, d, a],
            [b, d, a, c],
            [b, d, c, a],
            [c, a, b, d],
            [c, a, d, b],
            [c, b, a, d],
            [c, b, d, a],
            [c, d, a, b],
            [c, d, b, a],
            [d, a, b, c],
            [d, a, c, b],
            [d, b, a, c],
            [d, b, c, a],
            [d, c, a, b],
            [d, c, b, a],
        ]
    )
    weights = numpy.tile(w, 24)
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


def _s4_alt(dummy):
    return numpy.full((4, 1), 0.25)


def _s31_alt(a):
    b = 1 - 3 * a
    points = numpy.array([[a, a, a, b], [a, a, b, a], [a, b, a, a], [b, a, a, a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s22_alt(a):
    b = (1 - 2 * a) / 2
    points = numpy.array(
        [
            [a, a, b, b],
            [a, b, a, b],
            [b, a, a, b],
            [a, b, b, a],
            [b, a, b, a],
            [b, b, a, a],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s211_alt(data):
    a, b = data
    c = 1 - 2 * a - b
    points = numpy.array(
        [
            [a, a, b, c],
            [a, b, a, c],
            [b, a, a, c],
            [a, b, c, a],
            [b, a, c, a],
            [b, c, a, a],
            [a, a, c, b],
            [a, c, a, b],
            [c, a, a, b],
            [a, c, b, a],
            [c, a, b, a],
            [c, b, a, a],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s1111_alt(data):
    a, b, c = data
    d = 1 - a - b - c
    points = numpy.array(
        [
            [a, b, c, d],
            [a, b, d, c],
            [a, c, b, d],
            [a, c, d, b],
            [a, d, b, c],
            [a, d, c, b],
            [b, a, c, d],
            [b, a, d, c],
            [b, c, a, d],
            [b, c, d, a],
            [b, d, a, c],
            [b, d, c, a],
            [c, a, b, d],
            [c, a, d, b],
            [c, b, a, d],
            [c, b, d, a],
            [c, d, a, b],
            [c, d, b, a],
            [d, a, b, c],
            [d, a, c, b],
            [d, b, a, c],
            [d, b, c, a],
            [d, c, a, b],
            [d, c, b, a],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "s4": _s4_alt,
            "s31": _s31_alt,
            "s211": _s211_alt,
            "s22": _s22_alt,
            "s1111": _s1111_alt,
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


def _read(filepath, source):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return T3Scheme(name, weights, points, degree, source, tol)
