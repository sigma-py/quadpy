# -*- coding: utf-8 -*-
#
from .helpers import cartesian_to_spherical
from ..helpers import untangle, fsd


class HeoXu(object):
    '''
    Sangwoo Heo and Yuan Xu,
    Constructing Fully Symmetric Cubature Formulae for the Sphere,
    Mathematics of Computation,
    Vol. 70, No. 233 (Jan., 2001), pp. 269-279,
    <https://doi.org/10.1090/S0025-5718-00-01198-4>.

    Abstract:
    We construct symmetric cubature formulae of degrees in the 13-39 range for
    the surface measure on the unit sphere. We exploit a recently published
    correspondence between cubature formulae on the sphere and on the triangle.
    Specifically, a fully symmetric cubature formula for the surface measure on
    the unit sphere corresponds to a symmetric cubature formula for the
    triangle with weight function (u1u 2u 3)-1/2, where u1, u2, and u3 are
    homogeneous coordinates.
    '''
    def __init__(self, index):
        if index == '13':
            self.degree = 13
            data = [
                (0.013866592105, fsd(3, (1.0, 1))),
                (0.013050931863, fsd(3, (0.286640146767, 2), (0.914152532416, 1))),
                (0.013206423223, fsd(3, (0.659905001656, 2), (0.359236381200, 1))),
                (0.011942663555, fsd(3, (0.539490098706, 1), (0.841991943785, 1))),
                ]
        else:
            assert False

        self.points, self.weights = untangle(data)
        self.phi_theta = cartesian_to_spherical(self.points)
        return
