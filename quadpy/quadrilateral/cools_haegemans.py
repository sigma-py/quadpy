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
                numpy.full(8, 0.361130558151e-01),
                numpy.full(4, 0.535500902317e-01),
                numpy.full(4, 0.106828079664e-01),
                numpy.full(4, 0.113540990172)
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
                numpy.full(8, 0.348818790231e-01),
                numpy.full(8, 0.344998496602e-01),
                numpy.full(8, 0.987441946914e-02),
                numpy.full(8, 0.203490805188e-01),
                numpy.full(4, 0.475325029082e-01),
                numpy.full(4, 0.325703974952e-02),
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
            assert index == 3
            self.degree = 13
            self.weights = 4.0 * numpy.concatenate([
                numpy.full(8, 0.197386321888e-01),
                numpy.full(8, 0.484363166325e-01),
                numpy.full(8, 0.421281899422e-02),
                numpy.full(8, 0.287255968895e-01),
                numpy.full(4, 0.361061434781e-01),
                numpy.full(4, 0.116671271121e-01),
                ])
            self.points = numpy.concatenate([
                _s8(0.168234947696, 0.914794463441),
                _s8(0.252666976106, 0.591294378163),
                _s8(0.584047706043, 0.102139695463e+01),
                _s8(0.586713014973, 0.826081709475),
                _s4(0.178898689064),
                _s4(0.914197956909),
                ])
        # TODO There are three more schemes in the techincal report
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
