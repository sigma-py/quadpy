# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _fsd, _pm
from . import ewing
from . import hammer_stroud
from . import mustard_lyness_blatt
from . import phillips
from . import stroud57
from . import stroud66
from . import stroud68
from . import thacher
from . import tyler


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index):
        self.name = 'Stroud({})'.format(index)
        reference_volume = 2.0**n
        if index == 'Cn 1-1':
            # centroid formula
            self.degree = 1
            self.weights = numpy.array([reference_volume])
            self.points = numpy.array([
                numpy.full(n, 0.0)
                ])
        elif index == 'Cn 1-2':
            # product trapezoidal formula
            self.degree = 1
            self.weights = numpy.full(2**n, 1.0)
            self.points = _pm(n, 1.0)
        elif index == 'Cn 2-1':
            scheme = stroud57.Stroud57(n, 2)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 2-2':
            scheme = thacher.Thacher(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 3-1':
            scheme = stroud57.Stroud57(n, 3)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 3-2':
            self.degree = 3
            self.weights = numpy.full(2*n, reference_volume / (2*n))
            r = numpy.sqrt(n / 3.0)
            self.points = _fsd(n, r, 1)
        elif index == 'Cn 3-3':
            scheme = tyler.Tyler(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 3-4':
            # product Gauss formula
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(2**n, reference_volume / 2**n),
                ])
            r = numpy.sqrt(3.0) / 3.0
            self.points = _pm(n, r)
        elif index == 'Cn 3-5':
            scheme = ewing.Ewing(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 3-6':
            # product Simpson's formula
            self.degree = 3
            lst = n * [[1.0/3.0, 4.0/3.0, 1.0/3.0]]
            self.weights = numpy.product(
                numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n),
                axis=-1
                )
            lst = n * [[-1.0, 0.0, 1.0]]
            self.points = numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n)
        # elif index == 'Cn 5-1':
        # Cn 5-1 is not implemented because it's based on explicit values only
        # given for n=4,5,6.
        elif index == 'Cn 5-2':
            scheme = hammer_stroud.HammerStroud(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-3':
            scheme = stroud68.Stroud68(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-4':
            scheme = stroud66.Stroud66(n, 'a')
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-5':
            scheme = mustard_lyness_blatt.MustardLynessBlatt(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-6':
            scheme = stroud66.Stroud66(n, 'b')
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-7':
            scheme = stroud66.Stroud66(n, 'c')
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-8':
            scheme = stroud66.Stroud66(n, 'd')
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        elif index == 'Cn 5-9':
            # product Gauss formula
            self.degree = 5
            lst = n * [[5.0/9.0, 8.0/9.0, 5.0/9.0]]
            self.weights = numpy.product(
                numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n),
                axis=-1
                )
            sqrt35 = numpy.sqrt(3.0/5.0)
            lst = n * [[-sqrt35, 0.0, sqrt35]]
            self.points = numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n)
        else:
            assert index == 'Cn 7-1'
            scheme = phillips.Phillips(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights

        return
