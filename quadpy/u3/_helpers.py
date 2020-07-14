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
    return


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


def _a1():
    return numpy.array(
        [
            [+1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, +1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, +1.0],
            [0.0, 0.0, -1.0],
        ]
    ).T


def _a2():
    return numpy.array(
        [
            [+1.0, +1.0, 0.0],
            [+1.0, -1.0, 0.0],
            [-1.0, +1.0, 0.0],
            [-1.0, -1.0, 0.0],
            #
            [+1.0, 0.0, +1.0],
            [+1.0, 0.0, -1.0],
            [-1.0, 0.0, +1.0],
            [-1.0, 0.0, -1.0],
            #
            [0.0, +1.0, +1.0],
            [0.0, +1.0, -1.0],
            [0.0, -1.0, +1.0],
            [0.0, -1.0, -1.0],
        ]
    ).T / numpy.sqrt(2.0)


def _a3():
    return numpy.array(
        [
            [+1.0, +1.0, +1.0],
            [+1.0, +1.0, -1.0],
            [+1.0, -1.0, +1.0],
            [+1.0, -1.0, -1.0],
            [-1.0, +1.0, +1.0],
            [-1.0, +1.0, -1.0],
            [-1.0, -1.0, +1.0],
            [-1.0, -1.0, -1.0],
        ]
    ).T / numpy.sqrt(3.0)


def _pq0(alpha):
    # TODO merge with pq02
    a = numpy.sin(alpha * numpy.pi)
    b = numpy.cos(alpha * numpy.pi)
    zero = numpy.zeros_like(alpha)
    out = numpy.array(
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
    out = numpy.moveaxis(out, 0, 1)
    out = out.reshape(out.shape[0], -1)
    return out


def _pq02(a):
    b = numpy.sqrt(1 - a ** 2)
    zero = numpy.zeros_like(a)
    out = numpy.array(
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
    out = numpy.moveaxis(out, 0, 1)
    out = out.reshape(out.shape[0], -1)
    return out


def _llm(beta):
    # translate the point into cartesian coords; note that phi=pi/4.
    beta *= numpy.pi
    L = numpy.sin(beta) / numpy.sqrt(2)
    m = numpy.cos(beta)
    return _llm2(L, m)


def _llm2(L, m=None):
    if m is None:
        m = numpy.sqrt(1 - 2 * L ** 2)
    out = numpy.array(
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
    out = numpy.moveaxis(out, 0, 1)
    out = out.reshape(out.shape[0], -1)
    return out


def _rsw(azimuthal, polar):
    # translate the point into cartesian coords; note that phi=pi/4.
    azimuthal *= numpy.pi
    polar *= numpy.pi

    sin_polar = numpy.sin(polar)
    cos_polar = numpy.cos(polar)
    sin_azimuthal = numpy.sin(azimuthal)
    cos_azimuthal = numpy.cos(azimuthal)

    r = sin_polar * cos_azimuthal
    s = sin_polar * sin_azimuthal
    w = cos_polar

    return _rsw2(r, s, w)


def _rsw2(r, s, w=None):
    if w is None:
        w = numpy.sqrt(1 - r ** 2 - s ** 2)
    out = numpy.array(
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
    out = numpy.moveaxis(out, 0, 1)
    out = out.reshape(out.shape[0], -1)
    return out


def untangle2(data):
    points = []
    weights = []
    if "a1" in data:
        assert len(data["a1"]) == 1
        points.append(_a1())
        w = data["a1"][0]
        weights.append(numpy.full(6, w))

    if "a2" in data:
        assert len(data["a2"]) == 1
        points.append(_a2())
        w = data["a2"][0]
        weights.append(numpy.full(12, w))

    if "a3" in data:
        assert len(data["a3"]) == 1
        points.append(_a3())
        w = data["a3"][0]
        weights.append(numpy.full(8, w))

    if "llm" in data:
        llm = numpy.asarray(data["llm"]).T
        points.append(_llm(llm[1]))
        weights.append(numpy.tile(llm[0], 24))

    if "llm2" in data:
        llm = numpy.asarray(data["llm2"]).T
        points.append(_llm2(llm[1]))
        weights.append(numpy.tile(llm[0], 24))

    if "pq0" in data:
        pq0 = numpy.asarray(data["pq0"]).T
        points.append(_pq0(pq0[1]))
        weights.append(numpy.tile(pq0[0], 24))

    if "pq02" in data:
        pq0 = numpy.asarray(data["pq02"]).T
        points.append(_pq02(pq0[1]))
        weights.append(numpy.tile(pq0[0], 24))

    if "rsw" in data:
        rsw = numpy.asarray(data["rsw"]).T
        beta = rsw[1:]
        points.append(_rsw(*beta))
        weights.append(numpy.tile(rsw[0], 48))

    if "rsw2" in data:
        rsw = numpy.asarray(data["rsw2"]).T
        beta = rsw[1:]
        points.append(_rsw2(*beta))
        weights.append(numpy.tile(rsw[0], 48))

    if "plain" in data:
        dat = numpy.asarray(data["plain"])
        points.append(dat[:, :3].T)
        weights.append(dat[:, 3])

    points = numpy.ascontiguousarray(numpy.concatenate(points, axis=1))
    weights = numpy.concatenate(weights)
    return points, weights


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)


def _read(filepath, source, weight_factor=None):
    with open(filepath, "r") as f:
        content = json.load(f)

    degree = content.pop("degree")
    name = content.pop("name")
    tol = content.pop("test_tolerance")

    if tol > 1.0e-12:
        warnings.warn(f"The {name} scheme has low precision ({tol:.3e}).")

    points, weights = untangle2(content.pop("data"))
    theta_phi = cartesian_to_spherical(points)

    if weight_factor is not None:
        weights *= weight_factor

    return U3Scheme(name, weights, points, theta_phi, degree, source, tol=tol)
