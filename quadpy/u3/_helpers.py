import numpy
import sympy

from ..helpers import QuadratureScheme


class U3Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, theta_phi, degree, source):
        self.domain = "U3"
        self.name = name
        self.degree = degree
        self.source = source

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

        flt = numpy.vectorize(float)
        pts = flt(self.points)
        wgs = flt(self.weights)

        for p, w in zip(pts, wgs):
            # <https://en.wikipedia.org/wiki/Spherical_cap>
            w *= 4 * numpy.pi
            theta = numpy.arccos(1.0 - abs(w) / (2 * numpy.pi))
            color = "tab:blue" if w >= 0 else "tab:red"
            _plot_spherical_cap_mpl(ax, p, theta, color)

        ax.set_axis_off()

    def integrate(self, f, center, radius, dot=numpy.dot):
        """Quadrature where `f` is defined in Cartesian coordinates.
        """
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.array(f((rr + center).T))
        return area(radius) * dot(ff, self.weights)

    def integrate_spherical(self, f, dot=numpy.dot):
        """Quadrature where `f` is a function of the spherical coordinates theta_phi
        (polar, azimuthal, in this order).
        """
        ff = numpy.asarray(f(self.theta_phi.T))
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
    return numpy.stack([numpy.arccos(X[:, 2]), numpy.arctan2(X[:, 1], X[:, 0])], axis=1)


def _atan2_0(X):
    """Like sympy.atan2, but return 0 for x=y=0. Mathematically, the value is
    undefined, so sympy returns NaN, but for the sake of the coordinate
    conversion, its value doesn't matter. NaNs, however, produce NaNs down the
    line.
    """
    out = numpy.array([sympy.atan2(X[k, 1], X[k, 0]) for k in range(len(X))])
    out[out == sympy.nan] = 0
    return out


def cartesian_to_spherical_sympy(X):
    vacos = numpy.vectorize(sympy.acos)
    return numpy.stack([vacos(X[:, 2]), _atan2_0(X)], axis=1)


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
    )


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
    ) / numpy.sqrt(2.0)


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
    ) / numpy.sqrt(3.0)


def _pq0(alpha):
    a = numpy.sin(alpha * numpy.pi)
    b = numpy.cos(alpha * numpy.pi)
    zero = numpy.zeros_like(alpha)
    return numpy.array(
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


def _pq02(a):
    b = numpy.sqrt(1 - a ** 2)
    zero = numpy.zeros_like(a)
    return numpy.array(
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


def _llm(beta):
    # translate the point into cartesian coords; note that phi=pi/4.
    beta *= numpy.pi
    L = numpy.sin(beta) / numpy.sqrt(2)
    m = numpy.cos(beta)
    return numpy.array(
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


def _llm2(L):
    m = numpy.sqrt(1 - 2 * L ** 2)
    return numpy.array(
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

    return numpy.array(
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
        beta = numpy.array(data["llm"])[:, 1]
        out = _collapse0(numpy.moveaxis(_llm(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["llm"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "llm2" in data:
        beta = numpy.array(data["llm2"])[:, 1]
        out = _collapse0(numpy.moveaxis(_llm2(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["llm2"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "pq0" in data:
        beta = numpy.array(data["pq0"])[:, 1]
        out = _collapse0(numpy.moveaxis(_pq0(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["pq0"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "pq02" in data:
        beta = numpy.array(data["pq02"])[:, 1]
        out = _collapse0(numpy.moveaxis(_pq02(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["pq02"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "rsw" in data:
        beta = numpy.array(data["rsw"])[:, 1:].T
        out = _collapse0(numpy.moveaxis(_rsw(*beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["rsw"])[:, 0]
        weights.append(numpy.tile(w, 48))

    points = numpy.concatenate(points)
    weights = numpy.concatenate(weights)
    return points, weights


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)
