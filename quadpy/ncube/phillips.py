# -*- coding: utf-8 -*-
#
import warnings

import numpy

from .helpers import _fsd, _z


class Phillips(object):
    '''
    G.M. Phillips,
    Numerical integration over an N-dimensional rectangular region,
    Comput J (1967) 10 (3): 297-299,
    <https://doi.org/10.1093/comjnl/10.3.297>.

    Abstract:
    Gaussian-type formulae are derived for all values of N >= 2.
    '''
    def __init__(self, n):
        warnings.warn('The Phillips schemes are only single-precision.')

        self.name = 'Phillips'
        reference_volume = 2.0**n

        self.degree = 7

        if n == 2:
            r1 = 1.0
            r2 = 0.462910050
            s = 0.774596669
            B0 = -0.158024691
            B1 = 0.036363636
            B2 = 0.175982043
            B3 = 0.077160494
        elif n == 3:
            r1 = 1.0
            r2 = 0.597614305
            s = 0.632455532
            t = 1.0
            B0 = 0.542962963
            B1 = 0.032098765
            B2 = -0.193580247
            B3 = 0.115740741
            B4 = 0.004629630
        elif n == 4:
            r1 = 1.0
            r2 = 0.313391585
            s = 0.447213596
            t = 0.707106781
            B0 = -10.888215488
            B1 = 0.025082508
            B2 = 2.007240724
            B3 = -0.231481481
            B4 = 0.037037037
        elif n == 5:
            r1 = 1.0
            r2 = 0.572478028
            s = 0.774596669
            t = 0.774596669
            B0 = -0.987232247
            B1 = 0.042500000
            B2 = 0.190516778
            B3 = -0.051440329
            B4 = 0.021433471
        else:
            assert n == 6
            r1 = 1.0
            r2 = 0.684837715
            s = 0.774596669
            t = 0.774596669
            B0 = -1.620415307
            B1 = 0.053807107
            B2 = 0.350317580
            B3 = -0.094307270
            B4 = 0.021433471

        pts = [
            _z(n),
            _fsd(n, r1, 1),
            _fsd(n, r2, 1),
            _fsd(n, s, 2),
            ]
        wgts = [
            numpy.full(1, B0 * reference_volume),
            numpy.full(2*n, B1 * reference_volume),
            numpy.full(2*n, B2 * reference_volume),
            numpy.full(2*n*(n-1), B3 * reference_volume),
            ]

        if n > 2:
            pts.append(_fsd(n, t, 3))
            wgts.append(numpy.full(4*n*(n-1)*(n-2)//3, B4 * reference_volume))

        self.points = numpy.concatenate(pts)
        self.weights = numpy.concatenate(wgts)
        return
