# -*- coding: utf-8 -*-
#
from .helpers import untangle2


class MaeztuSainz(object):
    """
    J.I. Maeztu and E. Sainz de la Maza,
    An invariant quadrature rule of degree 11 for the tetrahedron,
    C. R. Acad. Sci. Paris 321 (1995), 1263-1267.
    """

    def __init__(self):
        self.name = "MaeztuSainz"
        # The article claims degree 11, but tests show only degree 1. :/
        # TODO find out what's going on
        self.degree = 1
        data = {
            "s4": [[-0.3229059250896649]],
            "s31": [
                [-0.3831136086645949, 0.3197881306061907],
                [+0.1259876832639002, 0.2745875432484354],
                [+0.7772656110490364e-2, 0.4902463231623282e-1],
                [+0.4475842042017354e-5, -0.58892050323316550e-1],
                [+0.3076630972851224e-1, 0.14369806508030763],
            ],
            "s22": [[+0.2230322290225118e-1, 0.43340593206769717]],
            "s211": [
                [+0.5167456484634155e-3, 0.50318342940324511, 0.60987466974805193e-1],
                [+0.1484538986489890, 0.29445616949492650, 0.37182110608410947],
                [+0.9330967352789100e-3, 0.0, 0.13838985309026736],
                [+0.9319130804165715e-2, 0.14550316358503807, 0.69267352508351802],
                [+0.1272850504266610e-1, 0.43854531792695007e-1, 0.27759599714708815],
            ],
        }

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
