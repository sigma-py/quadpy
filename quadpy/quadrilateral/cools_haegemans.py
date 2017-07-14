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
            self.degree = 9
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
        elif index == 2:
            self.degree = 13
            self.weights = 4.0 * numpy.concatenate([
                0.348818790231e-01 * numpy.ones(8),
                0.344998496602e-01 * numpy.ones(8),
                0.987441946914e-02 * numpy.ones(8),
                0.203490805188e-01 * numpy.ones(8),
                0.475325029082e-01 * numpy.ones(4),
                0.325703974952e-02 * numpy.ones(4),
                ])
            self.points = numpy.concatenate([
                _s8(0.266676738695e-01, 0.377724312590),
                _s8(0.235988332487, 0.793396171109),
                _s8(0.265486560241, 0.978761747825),
                _s8(0.702141598362, 0.913909457030),
                _s4(0.551473280570),
                _s4(0.968340720218),
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
