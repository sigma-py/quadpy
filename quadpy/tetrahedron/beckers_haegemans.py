# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s4, _s31, _s22, _s211


class BeckersHaegemans(object):
    '''
    M. Beckers and A. Haegemans,
    The construction of cubature formulae for the tetrahedron,
    Report TW 128, Dept. of Computer Science, K.U. Leuven, 1990,
    <https://lirias.kuleuven.be/handle/123456789/132648>.
    '''
    def __init__(self, index):
        self.name = 'BH({})'.format(index)
        if index == 8:
            self.degree = 8
            self.weights = 6 * numpy.concatenate([
                numpy.full(1, -0.020500188658639915),
                numpy.full(4, 0.014250305822866901),
                numpy.full(4, 1.9670333131339009e-3),
                numpy.full(4, 1.6983410909288737e-4),
                numpy.full(6, 4.5796838244672818e-3),
                numpy.full(12, 5.7044858086819185e-3),
                numpy.full(12, 2.1405191411620925e-3),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.20682993161067320),
                _s31(0.082103588310546723),
                _s31(5.7819505051979972e-3),
                _s22(0.44946725998110577),
                _s211(0.22906653611681113, 0.506227344977843697),
                _s211(0.036607749553197423, 0.19048604193463345),
                ])
        else:
            assert False

        self.points = bary[:, 1:]
        return
