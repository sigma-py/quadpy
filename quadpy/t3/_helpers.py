import json

import numpy as np
import sympy

from ..helpers import backend_to_function
from ..tn import TnScheme, get_vol, transform

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class T3Scheme(TnScheme):
    def __init__(
        self, name, symmetry_data, degree, source=None, tol=1.0e-14, weight_factor=None
    ):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data)
        if weight_factor is not None:
            weights *= weight_factor

        super().__init__(name, 3, weights, points, degree, source, tol)
        self.domain = "T3"

    def show(
        self,
        tet=np.array(
            [
                [+1, 0, -1.0 / np.sqrt(2.0)],
                [-1, 0, -1.0 / np.sqrt(2.0)],
                [0, +1, +1.0 / np.sqrt(2.0)],
                [0, -1, +1.0 / np.sqrt(2.0)],
            ]
        ),
        backend="vtk",
        render=True,
    ):
        edges = np.array([[tet[i], tet[j]] for i in range(4) for j in range(i)])
        edges = np.moveaxis(edges, 1, 2)
        backend_to_function[backend](
            transform(self.points.T, tet.T).T,
            self.weights,
            get_vol(tet),
            edges,
            render=render,
        )


def _s4(dummy):
    if dummy.dtype == sympy.Basic:
        return np.full((4, 1), sympy.Rational(1, 4))
    return np.full((4, 1), 0.25)


def _s31(a):
    b = 1 - 3 * a
    points = np.array([[a, a, a, b], [a, a, b, a], [a, b, a, a], [b, a, a, a]])
    points = np.moveaxis(points, 0, 1)
    return points


def _s22(a):
    if a.dtype in [sympy.Basic, int]:
        b = (1 - 2 * a) / sympy.S(2)
    else:
        b = (1 - 2 * a) / 2
    points = np.array(
        [
            [a, a, b, b],
            [a, b, a, b],
            [b, a, a, b],
            [a, b, b, a],
            [b, a, b, a],
            [b, b, a, a],
        ]
    )
    points = np.moveaxis(points, 0, 1)
    return points


def _s211(data):
    a, b = data
    c = 1 - 2 * a - b
    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def _s1111(data):
    a, b, c = data
    d = 1 - a - b - c
    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "s4": _s4,
            "s31": _s31,
            "s211": _s211,
            "s22": _s22,
            "s1111": _s1111,
            "plain": lambda vals: vals.reshape(4, 1, -1),
        }[key]
        pts = fun(np.asarray(points_raw))

        counts.append(pts.shape[1])
        pts = pts.reshape(pts.shape[0], -1)
        points.append(pts)

    points = np.ascontiguousarray(np.concatenate(points, axis=1))
    return points, counts


def expand_symmetries(data):
    # separate points and weights
    points_raw = {}
    weights_raw = []
    for key, values in data.items():
        values = np.asarray(values)
        weights_raw.append(values[0])
        points_raw[key] = values[1:]

    points, counts = expand_symmetries_points_only(points_raw)
    weights = np.concatenate(
        [np.tile(values, count) for count, values in zip(counts, weights_raw)]
    )
    return points, weights


def _read(filepath, source):
    with open(filepath) as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]
    data = content["data"]
    weight_factor = content["weight factor"] if "weight factor" in content else None
    return T3Scheme(name, data, degree, source, tol, weight_factor=weight_factor)
