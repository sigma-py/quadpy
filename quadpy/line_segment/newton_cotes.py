# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy


class NewtonCotesClosed(object):
    """
    Closed Newton-Cotes formulae.
    <https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas#Closed_Newton.E2.80.93Cotes_formulae>,
    <http://mathworld.wolfram.com/Newton-CotesFormulas.html>.
    """

    def __init__(self, index):
        self.points = numpy.linspace(-1.0, 1.0, index + 1)
        self.degree = index + 1 if index % 2 == 0 else index

        # Formula (26) from
        # <http://mathworld.wolfram.com/Newton-CotesFormulas.html>.
        # Note that Sympy carries out all operations in rationals, i.e.,
        # _exactly_. Only at the end, the rational is converted into a float.
        n = index
        self.weights = numpy.empty(n + 1)
        t = sympy.Symbol("t")
        for r in range(n + 1):
            # Compare with get_weights().
            f = sympy.prod([(t - i) for i in range(n + 1) if i != r])
            alpha = (
                2
                * (-1) ** (n - r)
                * sympy.integrate(f, (t, 0, n))
                / (math.factorial(r) * math.factorial(n - r))
                / index
            )
            self.weights[r] = alpha
        return


class NewtonCotesOpen(object):
    """
    Open Newton-Cotes formulae.
    <https://math.stackexchange.com/a/1959071/36678>
    """

    def __init__(self, index):
        self.points = numpy.linspace(-1.0, 1.0, index + 2)[1:-1]
        self.degree = index if (index + 1) % 2 == 0 else index - 1
        #
        n = index + 1
        self.weights = numpy.empty(n - 1)
        t = sympy.Symbol("t")
        for r in range(1, n):
            # Compare with get_weights().
            f = sympy.prod([(t - i) for i in range(1, n) if i != r])
            alpha = (
                2
                * (-1) ** (n - r + 1)
                * sympy.integrate(f, (t, 0, n))
                / (math.factorial(r - 1) * math.factorial(n - 1 - r))
                / n
            )
            self.weights[r - 1] = alpha
        return
