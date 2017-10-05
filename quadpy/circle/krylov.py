# -*- coding: utf-8 -*-
#
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
    def __init__(self, n):
        self.weights = numpy.full(n, 2*sympy.pi/n)
        self.points = numpy.column_stack([
            [sympy.cos(sympy.pi * sympy.Rational(2*k, n)) for k in range(n)],
            [sympy.sin(sympy.pi * sympy.Rational(2*k, n)) for k in range(n)],
            ])
        self.degree = n - 1
        return
