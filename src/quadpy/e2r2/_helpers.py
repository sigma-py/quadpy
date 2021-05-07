import json
import math

import numpy as np

from ..helpers import QuadratureScheme, expand_symmetries, plot_disks

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class E2r2Scheme(QuadratureScheme):
    def __init__(
        self, name, symmetry_data, degree, source, tol=1.0e-14, weight_factor=None
    ):
        self.domain = "E2r2"
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data, dim=2)
        assert points.shape[0] == 2
        if weight_factor is not None:
            weights *= weight_factor
        super().__init__(name, weights, points, degree, source, tol)

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        ax = plt.gca()
        plt.axis("equal")

        if not show_axes:
            ax.set_axis_off()

        I0 = 2 * math.pi
        plot_disks(plt, self.points.T, self.weights, I0)

    def integrate(self, f, dot=np.dot):
        flt = np.vectorize(float)
        ref_vol = math.pi
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))


def _read(filepath, source):
    with open(filepath) as f:
        content = json.load(f)

    degree = content["degree"]
    name = content["name"]
    tol = content["test_tolerance"]
    data = content["data"]

    weight_factor = content["weight factor"] if "weight factor" in content else None

    return E2r2Scheme(name, data, degree, source, tol, weight_factor=weight_factor)


def get_good_scheme(degree):
    if degree <= 15:
        return {
            0: schemes["stroud_4_1"],
            1: schemes["stroud_4_1"],
            2: schemes["stroud_4_1"],
            3: schemes["stroud_4_1"],
            4: schemes["stroud_4_1"],
            5: schemes["stroud_secrest_5"],
            6: schemes["stroud_secrest_6"],
            7: schemes["stroud_secrest_6"],
            8: schemes["haegemans_piessens_a"],
            9: schemes["haegemans_piessens_a"],
            10: schemes["rabinowitz_richter_2"],
            11: schemes["rabinowitz_richter_2"],
            12: schemes["rabinowitz_richter_5"],
            13: schemes["rabinowitz_richter_5"],
            14: schemes["rabinowitz_richter_5"],
            15: schemes["rabinowitz_richter_5"],
        }[degree]()

    return None
