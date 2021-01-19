import numpy as np
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
        points, weights = helpers.expand_symmetries(symmetry_data, dim=3)
        super().__init__(name, 3, weights, points, degree, source, tol)
        self.domain = "C3"

    def show(self, hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]), backend="vtk"):
        """Shows the quadrature points on a given hexahedron. The size of the balls
        around the points coincides with their weights."""
        edges = np.array(
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
        edges = np.moveaxis(edges, 1, 2)

        helpers.backend_to_function[backend](
            transform(self.points, hexa),
            self.weights,
            self.integrate(lambda x: np.ones(x.shape[1:]), hexa),
            edges,
        )


def get_good_scheme(degree):
    if degree <= 11:
        return {
            0: schemes["midpoint"](),
            1: schemes["midpoint"](),
            2: schemes["face_midpoint"](),
            3: schemes["face_midpoint"](),
            4: schemes["hammer_stroud_4_3"](),
            5: schemes["hammer_stroud_4_3"](),
            6: schemes["hammer_wymore"](sympy.Rational(27, 20)),
            7: schemes["hammer_wymore"](sympy.Rational(27, 20)),
            8: schemes["witherden_vincent_09"](),
            9: schemes["witherden_vincent_09"](),
            10: schemes["witherden_vincent_11"](),
            11: schemes["witherden_vincent_11"](),
        }[degree]

    return None
