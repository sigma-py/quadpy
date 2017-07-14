# -*- coding: utf-8 -*-
#
import numpy


class CoolsHaegemans(object):
    '''
    R. Cools, A. Haegemans,
    Construction of fully symmetric cubature formulae of degree 4k-3 for fully
    symmetric planar regions
    1985, Report TW 71, Dept. of Computer Science, KU Leuven,
    <https://lirias.kuleuven.be/bitstream/123456789/131870/1/TW71.pdf>.
    '''
    def __init__(self, index):
        self.name = 'CH(%d)' % index
        if index == 1:
            self.degree = 5
            self.weights = numpy.pi * numpy.concatenate([
                numpy.full(4, 0.233253175473),
                numpy.full(4, 0.167468245269e-01),
                ])
            self.points = numpy.concatenate([
                _s4(0.459700843381),
                _s40(0.125592606040e+01),
                ])
        elif index == 2:
            self.degree = 9
            self.weights = numpy.pi * numpy.concatenate([
                 numpy.full(8, 0.567209601536e-01),
                 numpy.full(4, 0.109948866164),
                 numpy.full(4, 0.261900192462e-01),
                 numpy.full(4, 0.419194282996e-03),
                ])
            self.points = numpy.concatenate([
                _s8(0.243244191752, 0.809458260086),
                _s4(0.302217386264),
                _s4(0.664341348594),
                _s40(0.134279080737e+01),
                ])
        else:
            assert index == 3
            self.degree = 9
            self.weights = numpy.pi * numpy.concatenate([
                 numpy.full(8, 0.123447696401e-01),
                 numpy.full(4, 0.932719633554e-01),
                 numpy.full(4, 0.589496783783e-01),
                 numpy.full(4, 0.730888189861e-01),
                ])
            self.points = numpy.concatenate([
                _s8(0.343855345294, 0.944778017142),
                _s4(0.277496500297),
                _s4(0.592355387396),
                _s40(0.778610819923),
                ])
        # TODO There are more schemes in the techincal report
        return


def _s8(a, b):
    return numpy.array([
        [+a, +b],
        [-a, +b],
        [+a, -b],
        [-a, -b],
        [+b, +a],
        [-b, +a],
        [+b, -a],
        [-b, -a],
        ])


def _s4(a):
    return numpy.array([
        [+a, +a],
        [-a, +a],
        [+a, -a],
        [-a, -a],
        ])


def _s40(a):
    return numpy.array([
        [+a, 0.0],
        [-a, 0.0],
        [0.0, +a],
        [0.0, -a],
        ])
