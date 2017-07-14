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
                0.233253175473 * numpy.ones(4),
                0.167468245269e-01 * numpy.ones(4),
                ])
            self.points = numpy.concatenate([
                _s4(0.459700843381),
                _s40(0.125592606040e+01),
                ])
        else:
            assert False
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
