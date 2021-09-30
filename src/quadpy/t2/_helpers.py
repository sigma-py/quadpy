from __future__ import annotations

import json

import numpy as np

from ..helpers import expand_symmetries, plot_disks
from ..tn import TnScheme, get_vol, transform

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class T2Scheme(TnScheme):
    def __init__(
        self,
        name: str,
        symmetry_data,
        degree: int,
        source=None,
        tol: float = 1.0e-14,
        comments: list[str] | None = None,
        weight_factor: float | None = None,
    ):
        self.symmetry_data = symmetry_data

        points, weights = expand_symmetries(symmetry_data, dim=2)
        if weight_factor is not None:
            weights *= weight_factor

        assert points.shape[0] == 3, f"{name}, {points.shape}"
        super().__init__(name, 2, weights, points, degree, source, tol, comments)
        self.domain = "T2"

    def plot(
        self,
        triangle=np.array([[-0.5, 0.0], [+0.5, 0.0], [0, 0.5 * (np.sqrt(3))]]),
        show_axes: bool = False,
    ):
        """Shows the quadrature points on a given triangle. The size of the circles
        around the points coincides with their weights.
        """
        from matplotlib import pyplot as plt

        plt.plot(triangle[:, 0], triangle[:, 1], "-k")
        plt.plot(
            [triangle[-1, 0], triangle[0, 0]], [triangle[-1, 1], triangle[0, 1]], "-k"
        )

        if not show_axes:
            plt.gca().set_axis_off()

        transformed_pts = transform(self.points, triangle.T).T

        vol = get_vol(triangle)
        plot_disks(plt, transformed_pts, self.weights, vol)

        plt.axis("equal")


def _read(filepath, source):
    with open(filepath) as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]
    d = content["data"]
    weight_factor = content["weight factor"] if "weight factor" in content else None
    return T2Scheme(name, d, degree, source, tol, weight_factor=weight_factor)


def _scheme_from_dict(content, source=None):
    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return T2Scheme(
        content["name"],
        weights,
        points,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )


def get_good_scheme(degree: int) -> T2Scheme | None:
    if degree <= 50:
        return {
            0: schemes["centroid"],
            1: schemes["centroid"],
            2: schemes["hammer_marlowe_stroud_3"],
            3: schemes["strang_fix_cowper_04"],
            4: schemes["williams_shunn_jameson_3"],
            5: schemes["hammer_marlowe_stroud_5"],
            6: schemes["strang_fix_cowper_09"],
            7: schemes["laursen_gellert_11"],
            8: schemes["lyness_jespersen_15"],
            9: schemes["lyness_jespersen_18"],
            10: schemes["witherden_vincent_10"],
            11: schemes["xiao_gimbutas_11"],
            12: schemes["xiao_gimbutas_12"],
            13: schemes["xiao_gimbutas_13"],
            14: schemes["xiao_gimbutas_14"],
            15: schemes["witherden_vincent_15"],
            16: schemes["xiao_gimbutas_16"],
            17: schemes["xiao_gimbutas_17"],
            18: schemes["witherden_vincent_18"],
            19: schemes["dunavant_19"],
            20: schemes["xiao_gimbutas_20"],
            21: schemes["xiao_gimbutas_21"],
            22: schemes["xiao_gimbutas_22"],
            23: schemes["xiao_gimbutas_23"],
            24: schemes["xiao_gimbutas_24"],
            25: schemes["xiao_gimbutas_25"],
            26: schemes["xiao_gimbutas_26"],
            27: schemes["xiao_gimbutas_27"],
            28: schemes["xiao_gimbutas_28"],
            29: schemes["xiao_gimbutas_29"],
            30: schemes["xiao_gimbutas_30"],
            31: schemes["xiao_gimbutas_31"],
            32: schemes["xiao_gimbutas_32"],
            33: schemes["xiao_gimbutas_33"],
            34: schemes["xiao_gimbutas_34"],
            35: schemes["xiao_gimbutas_35"],
            36: schemes["xiao_gimbutas_36"],
            37: schemes["xiao_gimbutas_37"],
            38: schemes["xiao_gimbutas_38"],
            39: schemes["xiao_gimbutas_39"],
            40: schemes["xiao_gimbutas_40"],
            41: schemes["xiao_gimbutas_41"],
            42: schemes["xiao_gimbutas_42"],
            43: schemes["xiao_gimbutas_43"],
            44: schemes["xiao_gimbutas_44"],
            45: schemes["xiao_gimbutas_45"],
            46: schemes["xiao_gimbutas_46"],
            47: schemes["xiao_gimbutas_47"],
            48: schemes["xiao_gimbutas_48"],
            49: schemes["xiao_gimbutas_49"],
            50: schemes["xiao_gimbutas_50"],
        }[degree]()

    # degree > 50
    return None
