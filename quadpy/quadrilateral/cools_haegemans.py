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
            # The article claims degree 9, but tests really only show degree 7.
            self.degree = 7
            self.weights = 4.0 * numpy.concatenate([
                0.361130558151e-01 * numpy.ones(8),
                0.535500902317e-01 * numpy.ones(4),
                0.106828079664e-01 * numpy.ones(4),
                0.113540990172 * numpy.ones(4)
                ])
            self.points = numpy.concatenate([
                _s8(0.344872025364, 0.918620441057),
                _s4(0.690880550486),
                _s4(0.939655258097),
                _s40(0.488926856974),
                ])
        else:
            assert False
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
