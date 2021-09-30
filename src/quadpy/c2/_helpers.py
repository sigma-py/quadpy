from __future__ import annotations

import json

import numpy as np

from ..cn import CnScheme
from ..cn import ncube_points as rectangle_points
from ..cn import transform
from ..helpers import expand_symmetries, plot_disks
from ..tn import get_vol

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class C2Scheme(CnScheme):
    def __init__(
        self,
        name: str,
        symmetry_data,
        degree: int,
        source=None,
        tol: float = 1.0e-14,
        comments: list[str] | None = None,
    ):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data, dim=2)
        super().__init__(name, 2, weights, points, degree, source, tol, comments)
        self.domain = "C2"

    def plot(self, quad=rectangle_points([0.0, 1.0], [0.0, 1.0]), show_axes=False):
        """Shows the quadrature points on a given quad. The area of the disks around the
        points coincides with their weights.
        """
        from matplotlib import pyplot as plt

        def plot_segment(a, b):
            plt.plot((a[0], b[0]), (a[1], b[1]), "-k")

        plot_segment(quad[0][0], quad[1][0])
        plot_segment(quad[1][0], quad[1][1])
        plot_segment(quad[1][1], quad[0][1])
        plot_segment(quad[0][1], quad[0][0])

        if not show_axes:
            plt.gca().set_axis_off()

        transformed_pts = transform(self.points, quad)

        # compute volume by splitting it in two triangles
        vol = get_vol(np.array([quad[0][0], quad[1][0], quad[0][1]])) + get_vol(
            np.array([quad[0][0], quad[0][1], quad[1][1]])
        )
        plot_disks(plt, transformed_pts, self.weights, vol)
        plt.axis("equal")
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 1.1)


def _read(filepath, source):
    with open(filepath) as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]

    data = content["data"]
    if "weight factor" in content:
        w = content["weight factor"]
        for value in data.values():
            value[0] = [val * w for val in value[0]]

    return C2Scheme(name, data, degree, source, tol)


def _scheme_from_dict(content, source=None):
    points, weights = expand_symmetries(content["data"])

    if "weight factor" in content:
        weights *= content["weight factor"]

    return C2Scheme(
        content["name"],
        weights,
        points,
        degree=content["degree"],
        source=source,
        tol=content["test_tolerance"],
        comments=content["comments"] if "comments" in content else None,
    )


def get_good_scheme(degree: int) -> C2Scheme | None:
    if degree <= 22:
        return {
            0: schemes["dunavant_00"],
            1: schemes["dunavant_00"],
            2: schemes["hammer_stroud_1_2"],
            3: schemes["hammer_stroud_1_2"],
            4: schemes["burnside"],
            5: schemes["burnside"],
            6: schemes["tyler_2"],
            7: schemes["tyler_2"],
            8: schemes["rabinowitz_richter_1"],
            9: schemes["rabinowitz_richter_1"],
            10: schemes["witherden_vincent_11"],
            11: schemes["witherden_vincent_11"],
            12: schemes["witherden_vincent_13"],
            13: schemes["witherden_vincent_13"],
            14: schemes["rabinowitz_richter_6"],
            15: schemes["rabinowitz_richter_6"],
            16: schemes["witherden_vincent_17"],
            17: schemes["witherden_vincent_17"],
            18: schemes["witherden_vincent_19"],
            19: schemes["witherden_vincent_19"],
            20: schemes["witherden_vincent_21"],
            21: schemes["witherden_vincent_21"],
        }[degree]()

    return None
