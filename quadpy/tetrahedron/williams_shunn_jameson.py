# -*- coding: utf-8 -*-
#
from .helpers import untangle2


class WilliamsShunnJameson(object):
    """
    D.M. Williams, L. Shunn, A. Jameson,
    Symmetric quadrature rules for simplexes based on sphere close packed
    lattice arrangements,
    Journal of Computational and Applied Mathematics,
    266 (2014) 18–38,
    <https://doi.org/10.1016/j.cam.2014.01.007>.

    Abstract:
    Sphere close packed (SCP) lattice arrangements of points are well-suited
    for formulating symmetric quadrature rules on simplexes, as they are
    symmetric under affine transformations of the simplex unto itself in 2D and
    3D. As a result, SCP lattice arrangements have been utilized to formulate
    symmetric quadrature rules with Np = 1, 4, 10, 20, 35, and 56 points on the
    3-simplex (Shunn and Ham, 2012). In what follows, the work on the 3-simplex
    is extended, and SCP lattices are employed to identify symmetric quadrature
    rules with Np = 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, and 66 points on the
    2-simplex and Np = 84 points on the 3-simplex. These rules are found to be
    capable of exactly integrating polynomials of up to degree 17 in 2D and up
    to degree 9 in 3D.
    """

    def __init__(self):
        self.degree = 9
        data = {
            "s31": [
                [0.002144935144316, 0.026878474414817],
                [0.020826641690769, 0.187140675803470],
                [0.023000681669286, 0.322111431830857],
            ],
            "s22": [
                [0.007210136064455, 0.473575835127937],
                [0.030798919159712, 0.352045262027356],
            ],
            "s211": [
                [0.004357844813864, 0.020953442220056, 0.225783205866940],
                [0.008593530677833, 0.096989733123466, 0.158462939666092],
                [0.004863063904912, 0.097608162890442, 0.011844417749498],
            ],
            "s1111": [
                [
                    0.015595140078259,
                    0.541184412800237,
                    0.133558160703568,
                    0.296501020543124,
                ]
            ],
        }
        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
