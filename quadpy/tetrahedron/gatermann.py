# -*- coding: utf-8 -*-
#
from .helpers import _s31, _s22

from ..helpers import untangle


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
        data = [
            (
                9.73033316198362119165356216965707e-06,
                _s31(0.656936552995394536166881327385593),
            ),
            (
                8.99031481668747219698547129902142e-03,
                _s31(0.0801424420792727848879183805550907),
            ),
            (
                2.17777476778781405656596945369837e-02,
                _s22(0.404475329343454044779549906725159),
            ),
        ]
        self.bary, self.weights = untangle(data)
        self.weights *= 6
        self.points = self.bary[:, 1:]
        return
