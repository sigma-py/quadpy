# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import untangle2


class Yu(object):
    """
    Yu Jinyun,
    Symmetyric Gaussian quadrature formulae for tetrahedronal regions,
    Computer Methods in Applied Mechanics and Engineering, 43 (1984) 349-353,
    <https://doi.org/10.1016/0045-7825(84)90072-0>.

    Abstract:
    Quadrature formulae of degrees 2 to 6 are presented for the numerical
    integration of a function over tetrahedronal regions. The formulae
    presented are of Gaussian type and fully symmetric with respect to the four
    vertices of the tetrahedron.
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "Yu({})".format(index)
        if index == 1:
            self.degree = 2
            data = {"s31": [[frac(1, 4), 0.138196601125015]]}
        elif index == 2:
            self.degree = 3
            data = {"s4": [[-frac(4, 5)]], "s31": [[frac(9, 20), frac(1, 6)]]}
        elif index == 3:
            self.degree = 4
            data = {
                "s31": [[0.5037379410012282e-01, 0.7611903264425430e-01]],
                "s211": [
                    [0.6654206863329239e-01, 0.4042339134672644, 0.1197005277978019]
                ],
            }
        elif index == 4:
            self.degree = 5
            data = {
                "s4": [[0.1884185567365411]],
                "s31": [[0.6703858372604275e-01, 0.8945436401412733e-01]],
                "s211": [
                    [0.4528559236327399e-01, 0.4214394310662522, 0.1325810999384657]
                ],
            }
        else:
            assert index == 5
            self.degree = 6
            data = {
                "s4": [[0.9040129046014750e-01]],
                "s31": [[0.1911983427899124e-01, 0.5742691731735682e-01]],
                "s211": [
                    [
                        0.4361493840666568e-01,
                        0.2312985436519147,
                        0.5135188412556341e-01,
                    ],
                    [
                        0.2581167596199161e-01,
                        0.4756909881472290e-01,
                        0.2967538129690260,
                    ],
                ],
            }

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
