# -*- coding: utf-8 -*-
#
import warnings

from .helpers import fs_r00, fs_rr0, pm_rrr

from ..helpers import untangle


class SarmaStroud(object):
    '''
    V.L.N. Sarma and A. H. Stroud,
    Eberlein Measure and Mechanical Quadrature Formulae. II. Numerical Results,
    Mathematics of Computation,
    Vol. 23, No. 108 (Oct., 1969), pp. 781-784,
    <https://dx.doi.org/10.2307/2004963>.
    '''
    def __init__(self):
        warnings.warn('SarmaStroud is only in single precision.')
        self.degree = 7

        data = [
            (0.3558180896e-1, fs_r00(0.9317380000)),
            (0.1247892770e-1, fs_rr0(0.9167441779)),
            (0.5286772991e-1, pm_rrr(0.4086003800)),
            (0.2672752182e-1, pm_rrr(0.7398529500)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
