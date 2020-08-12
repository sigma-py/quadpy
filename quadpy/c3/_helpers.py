import numpy

from .. import helpers
from ..cn import CnScheme
from ..cn import ncube_points as cube_points
from ..cn import transform


class C3Scheme(CnScheme):
    def __init__(self, name, weights, points, degree, source=None, tol=1.0e-14):
        assert points.shape[0] == 3
        super().__init__(name, 3, weights, points, degree, source, tol)
        self.domain = "C3"

    def show(self, hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]), backend="vtk"):
        """Shows the quadrature points on a given hexahedron. The size of the balls
        around the points coincides with their weights."""
        edges = numpy.array(
            [
                [hexa[0, 0, 0], hexa[1, 0, 0]],
                [hexa[1, 0, 0], hexa[1, 1, 0]],
                [hexa[1, 1, 0], hexa[0, 1, 0]],
                [hexa[0, 1, 0], hexa[0, 0, 0]],
                #
                [hexa[0, 0, 1], hexa[1, 0, 1]],
                [hexa[1, 0, 1], hexa[1, 1, 1]],
                [hexa[1, 1, 1], hexa[0, 1, 1]],
                [hexa[0, 1, 1], hexa[0, 0, 1]],
                #
                [hexa[0, 0, 0], hexa[0, 0, 1]],
                [hexa[1, 0, 0], hexa[1, 0, 1]],
                [hexa[1, 1, 0], hexa[1, 1, 1]],
                [hexa[0, 1, 0], hexa[0, 1, 1]],
            ]
        )
        edges = numpy.moveaxis(edges, 1, 2)

        helpers.backend_to_function[backend](
            transform(self.points, hexa),
            self.weights,
            self.integrate(lambda x: 1.0, hexa),
            edges,
        )


def z():
    return numpy.array([[0, 0, 0]])


def fs_r00(a):
    return numpy.array(
        [[+a, 0, 0], [0, +a, 0], [0, 0, +a], [-a, 0, 0], [0, -a, 0], [0, 0, -a]]
    )


def fs_rr0(a):
    return numpy.array(
        [
            [+a, +a, 0],
            [+a, 0, +a],
            [0, +a, +a],
            [+a, -a, 0],
            [+a, 0, -a],
            [0, +a, -a],
            [-a, +a, 0],
            [-a, 0, +a],
            [0, -a, +a],
            [-a, -a, 0],
            [-a, 0, -a],
            [0, -a, -a],
        ]
    )


def fs_rrs(a, b):
    return numpy.array(
        [
            [+a, +a, +b],
            [+a, +b, +a],
            [+b, +a, +a],
            [+a, -a, +b],
            [+a, +b, -a],
            [+b, +a, -a],
            [-a, +a, +b],
            [-a, +b, +a],
            [+b, -a, +a],
            [-a, -a, +b],
            [-a, +b, -a],
            [+b, -a, -a],
            [+a, +a, -b],
            [+a, -b, +a],
            [-b, +a, +a],
            [+a, -a, -b],
            [+a, -b, -a],
            [-b, +a, -a],
            [-a, +a, -b],
            [-a, -b, +a],
            [-b, -a, +a],
            [-a, -a, -b],
            [-a, -b, -a],
            [-b, -a, -a],
        ]
    )


def rss_pm(r, s):
    return numpy.array(
        [
            [+r, +s, +s],
            [+s, +r, +s],
            [+s, +s, +r],
            [-r, -s, -s],
            [-s, -r, -s],
            [-s, -s, -r],
        ]
    )


def pm_rrr(a):
    return numpy.array(
        [
            [+a, +a, +a],
            [-a, +a, +a],
            [+a, -a, +a],
            [-a, -a, +a],
            [+a, +a, -a],
            [-a, +a, -a],
            [+a, -a, -a],
            [-a, -a, -a],
        ]
    )


def _zero(data):
    return numpy.array([[0.0], [0.0], [0.0]])


def _symm_r00(r):
    zero = numpy.zeros_like(r)
    points = numpy.array([
        [+r, zero, zero], [-r, zero, zero],
        [zero, +r, zero], [zero, -r, zero],
        [zero, zero, +r], [zero, zero, -r]
    ])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_rr0(a):
    z = numpy.zeros_like(a)
    points = numpy.array([
        [+a, +a, z],
        [+a, z, +a],
        [z, +a, +a],
        [+a, -a, z],
        [+a, z, -a],
        [z, +a, -a],
        [-a, +a, z],
        [-a, z, +a],
        [z, -a, +a],
        [-a, -a, z],
        [-a, z, -a],
        [z, -a, -a],
    ])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_rrr(a):
    points = numpy.array([
        [+a, +a, +a],
        [-a, +a, +a],
        [+a, -a, +a],
        [-a, -a, +a],
        [+a, +a, -a],
        [-a, +a, -a],
        [+a, -a, -a],
        [-a, -a, -a],
    ])
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "zero": _zero,
            "symm_r00": _symm_r00,
            "symm_rr0": _symm_rr0,
            "symm_rrr": _symm_rrr,
            "plain": lambda vals: vals.reshape(3, 1, -1),
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
