# -*- coding: utf-8 -*-
#
import math

from .helpers import _s4, _s40, _s8
from ..helpers import untangle


class CoolsHaegemans(object):
    """
    R. Cools, A. Haegemans,
    Construction of fully symmetric cubature formulae of degree 4k-3 for fully
    symmetric planar regions
    1985, Report TW 71, Dept. of Computer Science, KU Leuven,
    <https://lirias.kuleuven.be/bitstream/123456789/131870/1/TW71.pdf>.
    """

    def __init__(self, index):
        self.name = "CH(%d)" % index
        if index == 1:
            self.degree = 5
            data = [
                (0.233253175473, _s4(0.459700843381)),
                (0.167468245269e-01, _s40(0.125592606040e01)),
            ]
        elif index == 2:
            self.degree = 9
            data = [
                (0.567209601536e-01, _s8(0.243244191752, 0.809458260086)),
                (0.109948866164e00, _s4(0.302217386264)),
                (0.261900192462e-01, _s4(0.664341348594)),
                (0.419194282996e-03, _s40(0.134279080737e01)),
            ]
        else:
            assert index == 3
            self.degree = 9
            data = [
                (0.123447696401e-01, _s8(0.343855345294, 0.944778017142)),
                (0.932719633554e-01, _s4(0.277496500297)),
                (0.589496783783e-01, _s4(0.592355387396)),
                (0.730888189861e-01, _s40(0.778610819923)),
            ]
        # TODO There are more schemes in the techincal report

        self.points, self.weights = untangle(data)
        self.weights *= math.pi
        return
