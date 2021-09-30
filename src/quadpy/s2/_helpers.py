from __future__ import annotations

import json
import warnings

import numpy as np

from .._exception import QuadpyError
from ..helpers import QuadratureScheme, expand_symmetries, plot_disks

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class S2Scheme(QuadratureScheme):
    def __init__(
        self,
        name: str,
        symmetry_data,
        degree: int,
        source=None,
        tol: float = 1.0e-14,
        weight_factor: float | None = None,
    ):
        self.domain = "S2"
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data, dim=2)
        assert points.shape[0] == 2
        assert points.shape[1] == weights.shape[0], f"{points.shape}, {weights.shape}"
        if weight_factor is not None:
            weights *= weight_factor
        super().__init__(name, weights, points, degree, source, tol)

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        ax = plt.gca()
        # change default range so that new disks will work
        plt.axis("equal")
        ax.set_xlim((-1.5, 1.5))
        ax.set_ylim((-1.5, 1.5))

        if not show_axes:
            ax.set_axis_off()

        disk1 = plt.Circle((0, 0), 1, color="k", fill=False)
        ax.add_patch(disk1)

        plot_disks(plt, self.points.T, self.weights, np.pi)

    def integrate(self, f, center, radius, dot=np.dot):
        center = np.array(center)
        rr = np.multiply.outer(radius, self.points.T)
        rr = np.swapaxes(rr, 0, -2)
        x = (rr + center).T
        fx = np.array(f(x))
        if fx.shape[-len(x.shape[1:]) :] != x.shape[1:]:
            string = ", ".join(str(val) for val in x.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {fx.shape}. Expected (..., {string})."
            )
        return np.pi * np.array(radius) ** 2 * dot(fx, self.weights)


def _read(filepath, source) -> S2Scheme:
    with open(filepath) as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    if "_ERR" in content:
        warnings.warn(content["_ERR"])

    weight_factor = content["weight factor"] if "weight factor" in content else None

    return S2Scheme(
        name, content["data"], degree, source, tol, weight_factor=weight_factor
    )


def _scheme_from_dict(content, source=None):
    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return S2Scheme(
        content["name"],
        weights,
        points,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )


def get_good_scheme(degree: int) -> S2Scheme | None:
    if degree <= 19:
        return {
            0: schemes["midpoint"],
            1: schemes["midpoint"],
            2: schemes["albrecht_collatz"],
            3: schemes["albrecht_collatz"],
            4: schemes["mysovskih_1"],
            5: schemes["radon"],
            6: schemes["kim_song_6"],
            7: schemes["kim_song_6"],
            8: schemes["luo_meng_2"],
            9: schemes["luo_meng_2"],
            10: schemes["mysovskih_2"],
            11: schemes["mysovskih_2"],
            12: schemes["cools_haegemans_13_1"],
            13: schemes["cools_haegemans_13_1"],
            14: schemes["kim_song_11"],
            15: schemes["kim_song_11"],
            16: schemes["cools_kim_1"],
            17: schemes["cools_kim_1"],
            18: schemes["kim_song_15"],
            19: schemes["kim_song_15"],
        }[degree]()

    return None
