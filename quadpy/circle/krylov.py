# -*- coding: utf-8 -*-
#
import numpy


class Krylov(object):
    '''
    Pages 73-74 in

    V.I. Krylov,
    Approximate Calculation of Integrals,
    Macmillan, New York, 1962.
    (Translated from 1st Russian ed., 1959, by A.H. Stroud.)
    <https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC&redir_esc=y>
    '''
    def __init__(self, n):
        self.weights = numpy.full(n, 2 * numpy.pi / n)
        self.points = numpy.column_stack([
            numpy.cos(2*numpy.pi * numpy.arange(n) / n),
            numpy.sin(2*numpy.pi * numpy.arange(n) / n),
            ])
        self.degree = n - 1
        return
