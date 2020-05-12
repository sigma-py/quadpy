import numpy
import sympy

from ..helpers import backend_to_function
from ..tn import TnScheme, get_vol, transform


class T3Scheme(TnScheme):
    def __init__(self, name, weights, points, degree, citation=None):
        self.name = name
        self.degree = degree
        self.citation = citation

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points
        return

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


def r(*data):
    w, r = numpy.array(data).T
    a = (1 - r) / 4
    b = 1 - 3 * a
    # like s31
    points = _stack_first_last([[a, a, a, b], [a, a, b, a], [a, b, a, a], [b, a, a, a]])
    weights = numpy.tile(w, 4)
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
