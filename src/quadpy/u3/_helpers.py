import json

import numpy as np
import orthopy
import sympy

from .._exception import QuadpyError
from ..helpers import QuadratureScheme

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class U3Scheme(QuadratureScheme):
    def __init__(
        self,
        name,
        symmetry_data,
        degree,
        source,
        tol=1.0e-14,
        comments=None,
    ):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data)

        if np.asarray(points).dtype == sympy.Basic:
            theta_phi = cartesian_to_spherical_sympy(points)
        else:
            theta_phi = cartesian_to_spherical(points)

        self.domain = "U3"
        super().__init__(name, weights, points, degree, source, tol, comments)

        theta_phi = np.asarray(theta_phi)
        if theta_phi.dtype == np.float64:
            self.theta_phi = theta_phi
        else:
            assert theta_phi.dtype in [np.dtype("O"), np.int_]
            self.theta_phi = theta_phi.astype(np.float64)
            self.theta_phi_symbolic = theta_phi

    def plot(self):
        from matplotlib import pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        ax = plt.axes(projection=Axes3D.name)
        # ax.set_aspect("equal")

        for p, w in zip(self.points.T, self.weights):
            # <https://en.wikipedia.org/wiki/Spherical_cap>
            w *= 4 * np.pi
            theta = np.arccos(1.0 - abs(w) / (2 * np.pi))
            color = "tab:blue" if w >= 0 else "tab:red"
            _plot_spherical_cap_mpl(ax, p, theta, color)

        ax.set_axis_off()
        return plt

    def integrate(self, f, center, radius, dot=np.dot):
        """Quadrature where `f` is defined in Cartesian coordinates."""
        center = np.asarray(center)
        rr = (self.points.T * radius + center).T
        frr = f(rr)
        if frr.shape[-len(rr.shape[1:]) :] != rr.shape[1:]:
            string = ", ".join(str(val) for val in rr.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {frr.shape}. Expected (..., {string})."
            )
        return area(radius) * dot(frr, self.weights)

    def integrate_spherical(self, f, dot=np.dot):
        """Quadrature where `f` is a function of the spherical coordinates theta_phi
        (polar, azimuthal, in this order).
        """
        ff = f(self.theta_phi)
        return area(1.0) * dot(ff, self.weights)

    def compute_residuals(self, level):
        evaluator = orthopy.u3.EvalCartesian(self.points, "quantum mechanic")

        max_res = []
        for k in range(level + 1):
            approximate = self.integrate(
                lambda x: next(evaluator), [0.0, 0.0, 0.0], 1.0
            )
            exact = evaluator.int_p0 if k == 0 else 0.0
            res = np.abs(approximate - exact)
            max_res += [np.max(res)]

        return np.array(max_res)


def area(radius):
    return 4 * np.pi * np.array(radius) ** 2


def _plot_spherical_cap_mpl(ax, b, opening_angle, color, elevation=1.01):
    r = elevation
    azimuthal = np.linspace(0, 2 * np.pi, 30)
    polar = np.linspace(0, opening_angle, 20)
    X = r * np.stack(
        [
            np.outer(np.cos(azimuthal), np.sin(polar)),
            np.outer(np.sin(azimuthal), np.sin(polar)),
            np.outer(np.ones(np.size(azimuthal)), np.cos(polar)),
        ],
        axis=-1,
    )

    # rotate X such that [0, 0, 1] gets rotated to `c`;
    # <https://math.stackexchange.com/a/476311/36678>.
    a = np.array([0.0, 0.0, 1.0])
    a_x_b = np.cross(a, b)
    a_dot_b = np.dot(a, b)
    if a_dot_b == -1.0:
        X_rot = -X
    else:
        X_rot = (
            X
            + np.cross(a_x_b, X)
            + np.cross(a_x_b, np.cross(a_x_b, X)) / (1.0 + a_dot_b)
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
    return np.array([np.arccos(X[2]), np.arctan2(X[1], X[0])])


def _atan2_0(X):
    """Like sympy.atan2, but return 0 for x=y=0. Mathematically, the value is
    undefined, so sympy returns NaN, but for the sake of the coordinate
    conversion, its value doesn't matter. NaNs, however, produce NaNs down the
    line.
    """
    out = np.array([sympy.atan2(X[1, k], X[0, k]) for k in range(X.shape[1])])
    out[out == sympy.nan] = 0
    return out


def cartesian_to_spherical_sympy(X):
    arccos = np.vectorize(sympy.acos)
    return np.array([arccos(X[2]), _atan2_0(X)])


def _a1(vals):
    symbolic = np.asarray(vals).dtype == sympy.Basic
    a = 1 if symbolic else 1.0
    points = np.array(
        [[+a, 0, 0], [-a, 0, 0], [0, +a, 0], [0, -a, 0], [0, 0, +a], [0, 0, -a]]
    ).T
    return points


def _a2(vals):
    symbolic = np.asarray(vals).dtype == sympy.Basic
    a = 1 / sympy.sqrt(2) if symbolic else 1 / np.sqrt(2)
    points = np.array(
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
    symbolic = np.asarray(vals).dtype == sympy.Basic
    a = 1 / sympy.sqrt(3) if symbolic else 1 / np.sqrt(3)
    points = np.array(
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
    return _pq02([np.sin(vals[0] * np.pi), np.cos(vals[0] * np.pi)])


def _pq02(vals):
    if len(vals) == 1:
        a = vals[0]
        b = np.sqrt(1 - a ** 2)
    else:
        assert len(vals) == 2
        a, b = vals

    if isinstance(a, sympy.Basic):
        zero = 0
    else:
        zero = np.zeros_like(a)

    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def _rs0(vals):
    if len(vals) == 1:
        a = vals[0]
        b = np.sqrt(1 - a ** 2)
    else:
        assert len(vals) == 2
        a, b = vals

    if isinstance(a, sympy.Basic):
        zero = 0
    else:
        zero = np.zeros_like(a)

    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def _llm(vals):
    # translate the point into cartesian coords; note that phi=pi/4.
    beta = vals[0] * np.pi
    L = np.sin(beta) / np.sqrt(2)
    m = np.cos(beta)
    return _llm2([L, m])


def _llm2(vals):
    if len(vals) == 1:
        L = vals[0]
        m = np.sqrt(1 - 2 * L ** 2)
    else:
        assert len(vals) == 2
        L, m = vals

    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def _rsw(vals):
    # translate the point into cartesian coords; note that phi=pi/4.
    phi_theta = vals * np.pi

    sin_phi, sin_theta = np.sin(phi_theta)
    cos_phi, cos_theta = np.cos(phi_theta)

    return _rsw2([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta])


def _rsw2(vals):
    if len(vals) == 2:
        r, s = vals
        w = np.sqrt(1 - r ** 2 - s ** 2)
    else:
        assert len(vals) == 3
        r, s, w = vals

    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def _rst(vals):
    if len(vals) == 2:
        r, s = vals
        w = np.sqrt(1 - r ** 2 - s ** 2)
    else:
        assert len(vals) == 3
        r, s, w = vals

    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
    return points


def _rst_weird(vals):
    if len(vals) == 2:
        r, s = vals
        t = np.sqrt(1 - r ** 2 - s ** 2)
    else:
        assert len(vals) == 3
        r, s, t = vals

    points = np.array(
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
    points = np.moveaxis(points, 0, 1)
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
        weights_raw.append(values[0])
        points_raw[key] = values[1:]

    points, counts = expand_symmetries_points_only(points_raw)
    weights = np.concatenate(
        [np.tile(values, count) for count, values in zip(counts, weights_raw)]
    )
    return points, weights


def _scheme_from_dict(content, source=None):
    data = content["data"]
    if "weight factor" in content:
        w = content["weight factor"]
        for val in data.values():
            val[0] = [v * w for v in val[0]]

    return U3Scheme(
        content["name"],
        data,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )


def _read(filepath, source):
    with open(filepath) as f:
        content = json.load(f)
    return _scheme_from_dict(content, source)


def get_good_scheme(degree):
    if degree <= 47:
        return {
            0: schemes["lebedev_003a"],
            1: schemes["lebedev_003a"],
            2: schemes["lebedev_003a"],
            3: schemes["lebedev_003a"],
            4: schemes["albrecht_collatz_1"],
            5: schemes["albrecht_collatz_1"],
            6: schemes["albrecht_collatz_5"],
            7: schemes["albrecht_collatz_5"],
            8: schemes["mclaren_05"],
            9: schemes["mclaren_05"],
            10: schemes["mclaren_08"],
            11: schemes["mclaren_08"],
            12: schemes["lebedev_015"],
            13: schemes["lebedev_015"],
            14: schemes["lebedev_015"],
            15: schemes["lebedev_015"],
            16: schemes["lebedev_017"],
            17: schemes["lebedev_017"],
            18: schemes["lebedev_019"],
            19: schemes["lebedev_019"],
            20: schemes["lebedev_021"],
            21: schemes["lebedev_021"],
            22: schemes["lebedev_023"],
            23: schemes["lebedev_023"],
            24: schemes["heo_xu_25b"],
            25: schemes["heo_xu_25b"],
            26: schemes["lebedev_029"],
            27: schemes["lebedev_029"],
            28: schemes["lebedev_029"],
            29: schemes["lebedev_029"],
            30: schemes["lebedev_031"],
            31: schemes["lebedev_031"],
            32: schemes["lebedev_035"],
            33: schemes["lebedev_035"],
            34: schemes["lebedev_035"],
            35: schemes["lebedev_035"],
            36: schemes["lebedev_041"],
            37: schemes["lebedev_041"],
            38: schemes["lebedev_041"],
            39: schemes["lebedev_041"],
            40: schemes["lebedev_041"],
            41: schemes["lebedev_041"],
            42: schemes["lebedev_047"],
            43: schemes["lebedev_047"],
            44: schemes["lebedev_047"],
            45: schemes["lebedev_047"],
            46: schemes["lebedev_047"],
            47: schemes["lebedev_047"],
        }[degree]()

    return None
