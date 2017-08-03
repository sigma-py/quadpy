# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _z, _symm_s_t


class Maxwell(object):
    '''
    J.C. Maxwell,
    On Approximate Multiple Integration between Limits by Summation.
    In W. Niven (Ed.), The Scientific Papers of James Clerk Maxwell,
    Cambridge Library Collection - Physical Sciences, pp. 604-611.
    Cambridge: Cambridge University Press.
    <https://doi.org/10.1017/CBO9780511710377.061>.
    '''
    def __init__(self):
        self.name = 'Maxwell'
        self.degree = 7

        r = numpy.sqrt(12.0 / 35.0)
        s = numpy.sqrt((93.0 + 3.0*numpy.sqrt(186.0)) / 155.0)
        t = numpy.sqrt((93.0 - 3.0*numpy.sqrt(186.0)) / 155.0)

        self.weights = numpy.concatenate([
            numpy.full(1, 1.0/81.0),
            numpy.full(4, 49.0/324.0),
            # typo in Stroud: 648 vs 649
            numpy.full(8, 31.0/648.0),
            ])
        self.points = numpy.concatenate([
            _z(),
            _symm_r_0(r),
            _symm_s_t(s, t)
            ])

        self.weights *= 4.0
        return
