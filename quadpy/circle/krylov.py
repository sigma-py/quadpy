# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy


class Krylov(object):
    '''
    Pages 73-74 in

    V.I. Krylov,
    Approximate Calculation of Integrals,
    Macmillan, New York, 1962.
    (Translated from 1st Russian ed., 1959, by A.H. Stroud.)
    <https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC&redir_esc=y>
    '''
    def __init__(self, n, symbolic=False):
        cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
        sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
        pi = sympy.pi if symbolic else numpy.pi

        self.weights = numpy.full(n, 2*pi/n)
        alpha = 2*numpy.arange(n) * pi / n
        self.points = numpy.column_stack([cos(alpha), sin(alpha)])
        self.degree = n - 1
        return
