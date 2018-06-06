# -*- coding: utf-8 -*-
#
from .helpers import untangle2


class Gatermann(object):
    """
    Karin Gatermann,
    Linear Representations of Finite Groups and The Ideal Theoretical
    Construction of G-Invariant Cubature Formulas,
    Numerical Integration pp 25-35,
    Part of the NATO ASI Series book series (ASIC, volume 357).
    """

    def __init__(self):
        self.name = "Gatermann"
        self.degree = 5
        data = {
            "s31": [
                [
                    9.73033316198362119165356216965707e-06,
                    0.656936552995394536166881327385593,
                ],
                [
                    8.99031481668747219698547129902142e-03,
                    0.0801424420792727848879183805550907,
                ],
            ],
            "s22": [
                [
                    2.17777476778781405656596945369837e-02,
                    0.404475329343454044779549906725159,
                ]
            ],
        }
        self.bary, self.weights = untangle2(data)
        self.weights *= 6
        self.points = self.bary[:, 1:]
        return
