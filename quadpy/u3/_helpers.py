import json
import warnings

import numpy
import sympy

from ..helpers import QuadratureScheme


class U3Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, theta_phi, degree, source, tol=1.0e-14):
        self.domain = "U3"
        self.name = name
        self.degree = degree
        self.source = source
        self.test_tolerance = tol

        weights = numpy.asarray(weights)
        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        points = numpy.asarray(points)
        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points

        assert weights.shape[0] == points.shape[1], (
            f"Shape mismatch for {name}: "
            f"weights.shape = {weights.shape}, points.shape = {points.shape}"
        )

        theta_phi = numpy.asarray(theta_phi)
        if theta_phi.dtype == numpy.float64:
            self.theta_phi = theta_phi
        else:
            assert theta_phi.dtype in [numpy.dtype("O"), numpy.int_]
            self.theta_phi = theta_phi.astype(numpy.float64)
            self.theta_phi_symbolic = theta_phi

    def plot(self):
        from matplotlib import pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.gca(projection=Axes3D.name)
        # ax.set_aspect("equal")

        for p, w in zip(self.points.T, self.weights):
            # <https://en.wikipedia.org/wiki/Spherical_cap>
            w *= 4 * numpy.pi
            theta = numpy.arccos(1.0 - abs(w) / (2 * numpy.pi))
            color = "tab:blue" if w >= 0 else "tab:red"
            _plot_spherical_cap_mpl(ax, p, theta, color)

        ax.set_axis_off()

    def integrate(self, f, center, radius, dot=numpy.dot):
        """Quadrature where `f` is defined in Cartesian coordinates.
        """
        center = numpy.asarray(center)
        rr = (self.points.T * radius + center).T
        return area(radius) * dot(f(rr), self.weights)

    def integrate_spherical(self, f, dot=numpy.dot):
        """Quadrature where `f` is a function of the spherical coordinates theta_phi
        (polar, azimuthal, in this order).
        """
        ff = f(self.theta_phi)
        return area(1.0) * dot(ff, self.weights)


def area(radius):
    return 4 * numpy.pi * numpy.array(radius) ** 2


def _plot_spherical_cap_mpl(ax, b, opening_angle, color, elevation=1.01):
    r = elevation
    azimuthal = numpy.linspace(0, 2 * numpy.pi, 30)
    polar = numpy.linspace(0, opening_angle, 20)
    X = r * numpy.stack(
        [
            numpy.outer(numpy.cos(azimuthal), numpy.sin(polar)),
            numpy.outer(numpy.sin(azimuthal), numpy.sin(polar)),
            numpy.outer(numpy.ones(numpy.size(azimuthal)), numpy.cos(polar)),
        ],
        axis=-1,
    )

    # rotate X such that [0, 0, 1] gets rotated to `c`;
    # <https://math.stackexchange.com/a/476311/36678>.
    a = numpy.array([0.0, 0.0, 1.0])
    a_x_b = numpy.cross(a, b)
    a_dot_b = numpy.dot(a, b)
    if a_dot_b == -1.0:
        X_rot = -X
    else:
        X_rot = (
            X
            + numpy.cross(a_x_b, X)
            + numpy.cross(a_x_b, numpy.cross(a_x_b, X)) / (1.0 + a_dot_b)
        )

    ax.plot_surface(
        X_rot[..., 0],
        X_rot[..., 1],
        X_rot[..., 2],
        rstride=3,
        cstride=3,
        color=color,
        alpha=0.5,
        linewidth=0,
    )


def cartesian_to_spherical(X):
    return numpy.array([numpy.arccos(X[2]), numpy.arctan2(X[1], X[0])])


def _atan2_0(X):
    """Like sympy.atan2, but return 0 for x=y=0. Mathematically, the value is
    undefined, so sympy returns NaN, but for the sake of the coordinate
    conversion, its value doesn't matter. NaNs, however, produce NaNs down the
    line.
    """
    out = numpy.array([sympy.atan2(X[1, k], X[0, k]) for k in range(X.shape[1])])
    out[out == sympy.nan] = 0
    return out


def cartesian_to_spherical_sympy(X):
    arccos = numpy.vectorize(sympy.acos)
    return numpy.array([arccos(X[2]), _atan2_0(X)])


def _a1(vals):
    symbolic = numpy.asarray(vals).dtype == sympy.Basic
    a = 1 if symbolic else 1.0
    points = numpy.array(
        [[+a, 0, 0], [-a, 0, 0], [0, +a, 0], [0, -a, 0], [0, 0, +a], [0, 0, -a]]
    ).T
    return points


def _a2(vals):
    symbolic = numpy.asarray(vals).dtype == sympy.Basic
    a = 1 / sympy.sqrt(2) if symbolic else 1 / numpy.sqrt(2)
    points = numpy.array(
        [
            [+a, +a, 0],
            [+a, -a, 0],
            [-a, +a, 0],
            [-a, -a, 0],
            #
            [+a, 0, +a],
            [+a, 0, -a],
            [-a, 0, +a],
            [-a, 0, -a],
            #
            [0, +a, +a],
            [0, +a, -a],
            [0, -a, +a],
            [0, -a, -a],
        ]
    ).T
    return points


def _a3(vals):
    symbolic = numpy.asarray(vals).dtype == sympy.Basic
    a = 1 / sympy.sqrt(3) if symbolic else 1 / numpy.sqrt(3)
    points = numpy.array(
        [
            [+a, +a, +a],
            [+a, +a, -a],
            [+a, -a, +a],
            [+a, -a, -a],
            [-a, +a, +a],
            [-a, +a, -a],
            [-a, -a, +a],
            [-a, -a, -a],
        ]
    ).T
    return points


def _pq0(vals):
    return _pq02([numpy.sin(vals[0] * numpy.pi), numpy.cos(vals[0] * numpy.pi)])


def _pq02(vals):
    if len(vals) == 1:
        a = vals[0]
        b = numpy.sqrt(1 - a ** 2)
    else:
        assert len(vals) == 2
        a, b = vals

    if isinstance(a, sympy.Basic):
        zero = 0
    else:
        zero = numpy.zeros_like(a)

    points = numpy.array(
        [
            [+a, +b, zero],
            [-a, +b, zero],
            [-a, -b, zero],
            [+a, -b, zero],
            #
            [+b, +a, zero],
            [-b, +a, zero],
            [-b, -a, zero],
            [+b, -a, zero],
            #
            [+a, zero, +b],
            [-a, zero, +b],
            [-a, zero, -b],
            [+a, zero, -b],
            #
            [+b, zero, +a],
            [-b, zero, +a],
            [-b, zero, -a],
            [+b, zero, -a],
            #
            [zero, +a, +b],
            [zero, -a, +b],
            [zero, -a, -b],
            [zero, +a, -b],
            #
            [zero, +b, +a],
            [zero, -b, +a],
            [zero, -b, -a],
            [zero, +b, -a],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _rs0(vals):
    if len(vals) == 1:
        a = vals[0]
        b = numpy.sqrt(1 - a ** 2)
    else:
        assert len(vals) == 2
        a, b = vals

    if isinstance(a, sympy.Basic):
        zero = 0
    else:
        zero = numpy.zeros_like(a)

    points = numpy.array(
        [
            [+a, +b, zero],
            [-a, +b, zero],
            [-a, -b, zero],
            [+a, -b, zero],
            #
            [+b, zero, +a],
            [-b, zero, +a],
            [-b, zero, -a],
            [+b, zero, -a],
            #
            [zero, +a, +b],
            [zero, -a, +b],
            [zero, -a, -b],
            [zero, +a, -b],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _llm(vals):
    # translate the point into cartesian coords; note that phi=pi/4.
    beta = vals[0] * numpy.pi
    L = numpy.sin(beta) / numpy.sqrt(2)
    m = numpy.cos(beta)
    return _llm2([L, m])


def _llm2(vals):
    if len(vals) == 1:
        L = vals[0]
        m = numpy.sqrt(1 - 2 * L ** 2)
    else:
        assert len(vals) == 2
        L, m = vals

    points = numpy.array(
        [
            [+L, +L, +m],
            [-L, +L, +m],
            [+L, -L, +m],
            [-L, -L, +m],
            [+L, +L, -m],
            [-L, +L, -m],
            [+L, -L, -m],
            [-L, -L, -m],
            #
            [+L, +m, +L],
            [-L, +m, +L],
            [+L, +m, -L],
            [-L, +m, -L],
            [+L, -m, +L],
            [-L, -m, +L],
            [+L, -m, -L],
            [-L, -m, -L],
            #
            [+m, +L, +L],
            [+m, -L, +L],
            [+m, +L, -L],
            [+m, -L, -L],
            [-m, +L, +L],
            [-m, -L, +L],
            [-m, +L, -L],
            [-m, -L, -L],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _rsw(vals):
    # translate the point into cartesian coords; note that phi=pi/4.
    phi_theta = vals * numpy.pi

    sin_phi, sin_theta = numpy.sin(phi_theta)
    cos_phi, cos_theta = numpy.cos(phi_theta)

    return _rsw2([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta])


def _rsw2(vals):
    if len(vals) == 2:
        r, s = vals
        w = numpy.sqrt(1 - r ** 2 - s ** 2)
    else:
        assert len(vals) == 3
        r, s, w = vals

    points = numpy.array(
        [
            [+r, +s, +w],
            [+w, +r, +s],
            [+s, +w, +r],
            [+s, +r, +w],
            [+w, +s, +r],
            [+r, +w, +s],
            #
            [-r, +s, +w],
            [+w, -r, +s],
            [+s, +w, -r],
            [+s, -r, +w],
            [+w, +s, -r],
            [-r, +w, +s],
            #
            [+r, -s, +w],
            [+w, +r, -s],
            [-s, +w, +r],
            [-s, +r, +w],
            [+w, -s, +r],
            [+r, +w, -s],
            #
            [+r, +s, -w],
            [-w, +r, +s],
            [+s, -w, +r],
            [+s, +r, -w],
            [-w, +s, +r],
            [+r, -w, +s],
            #
            [-r, -s, +w],
            [+w, -r, -s],
            [-s, +w, -r],
            [-s, -r, +w],
            [+w, -s, -r],
            [-r, +w, -s],
            #
            [-r, +s, -w],
            [-w, -r, +s],
            [+s, -w, -r],
            [+s, -r, -w],
            [-w, +s, -r],
            [-r, -w, +s],
            #
            [+r, -s, -w],
            [-w, +r, -s],
            [-s, -w, +r],
            [-s, +r, -w],
            [-w, -s, +r],
            [+r, -w, -s],
            #
            [-r, -s, -w],
            [-w, -r, -s],
            [-s, -w, -r],
            [-s, -r, -w],
            [-w, -s, -r],
            [-r, -w, -s],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _rst(vals):
    if len(vals) == 2:
        r, s = vals
        w = numpy.sqrt(1 - r ** 2 - s ** 2)
    else:
        assert len(vals) == 3
        r, s, w = vals

    points = numpy.array(
        [
            [+r, +s, +w],
            [+w, +r, +s],
            [+s, +w, +r],
            #
            [-r, +s, +w],
            [+w, -r, +s],
            [+s, +w, -r],
            #
            [+r, -s, +w],
            [+w, +r, -s],
            [-s, +w, +r],
            #
            [+r, +s, -w],
            [-w, +r, +s],
            [+s, -w, +r],
            #
            [-r, -s, +w],
            [+w, -r, -s],
            [-s, +w, -r],
            #
            [-r, +s, -w],
            [-w, -r, +s],
            [+s, -w, -r],
            #
            [+r, -s, -w],
            [-w, +r, -s],
            [-s, -w, +r],
            #
            [-r, -s, -w],
            [-w, -r, -s],
            [-s, -w, -r],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _rst_weird(vals):
    if len(vals) == 2:
        r, s = vals
        t = numpy.sqrt(1 - r ** 2 - s ** 2)
    else:
        assert len(vals) == 3
        r, s, t = vals

    points = numpy.array(
        [
            [+r, +s, +t],
            [-r, +t, +s],
            [+s, +t, +r],
            [-s, +r, +t],
            [+t, +r, +s],
            [-t, +s, +r],
            #
            [+r, -s, -t],
            [-r, -t, -s],
            [+s, -t, -r],
            [-s, -r, -t],
            [+t, -r, -s],
            [-t, -s, -r],
            #
            [+r, +t, -s],
            [-r, +s, -t],
            [+s, +r, -t],
            [-s, +t, -r],
            [+t, +s, -r],
            [-t, +r, -s],
            #
            [+r, -t, +s],
            [-r, -s, +t],
            [+s, -r, +t],
            [-s, -t, +r],
            [+t, -s, +r],
            [-t, -r, +s],
        ]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "a1": _a1,
            "a2": _a2,
            "a3": _a3,
            "llm": _llm,
            "llm2": _llm2,
            "pq0": _pq0,
            "pq02": _pq02,
            "rs0": _rs0,
            "rsw": _rsw,
            "rsw2": _rsw2,
            "rst": _rst,
            "rst_weird": _rst_weird,
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


def _read(filepath, source, weight_factor=None):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content.pop("degree")
    name = content.pop("name")
    tol = content.pop("test_tolerance")

    if tol > 1.0e-12:
        warnings.warn(f"The {name} scheme has low precision ({tol:.3e}).")

    points, weights = expand_symmetries(content.pop("data"))
    theta_phi = cartesian_to_spherical(points)

    if weight_factor is not None:
        weights *= weight_factor

    return U3Scheme(name, weights, points, theta_phi, degree, source, tol=tol)
