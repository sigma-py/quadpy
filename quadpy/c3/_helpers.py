import numpy
import sympy

from .. import helpers
from ..cn import CnScheme
from ..cn import ncube_points as cube_points
from ..cn import transform

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class C3Scheme(CnScheme):
    def __init__(self, name, symmetry_data, degree, source=None, tol=1.0e-14):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data)
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


def _zero(data):
    return numpy.array([[0.0], [0.0], [0.0]])


def _symm_r00(r):
    zero = numpy.zeros_like(r)
    points = numpy.array(
        [
            [+r, zero, zero],
            [-r, zero, zero],
            [zero, +r, zero],
            [zero, -r, zero],
            [zero, zero, +r],
            [zero, zero, -r],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_rr0(a):
    z = numpy.zeros_like(a)
    points = numpy.array(
        [
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
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_rrr(a):
    points = numpy.array(
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
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_rrs(data):
    a, b = data
    points = numpy.array(
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
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_rss_pm(data):
    r, s = data
    points = numpy.array(
        [
            [+r, +s, +s],
            [+s, +r, +s],
            [+s, +s, +r],
            [-r, -s, -s],
            [-s, -r, -s],
            [-s, -s, -r],
        ]
    )
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
            "symm_rrs": _symm_rrs,
            "symm_rss_pm": _symm_rss_pm,
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


def get_good_scheme(degree):
    if degree <= 7:
        return {
            0: schemes["midpoint"](),
            1: schemes["midpoint"](),
            2: schemes["face_midpoint"](),
            3: schemes["face_midpoint"](),
            4: schemes["hammer_stroud_4_3"](),
            5: schemes["hammer_stroud_4_3"](),
            6: schemes["hammer_wymore"](sympy.Rational(27, 20)),
            7: schemes["hammer_wymore"](sympy.Rational(27, 20)),
        }[degree]

    return None
