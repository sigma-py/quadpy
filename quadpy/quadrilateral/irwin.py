# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import concat, symm_s, symm_s_t


class Irwin(object):
    """
    Joseph Oscar Irwin,
    On quadrature and cubature,
    Cambridge University Press, 1923,
    <https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ&redir_esc=y>
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "Irwin({})".format(index)
        if index == 1:
            self.degree = 3
            self.weights, self.points = concat(
                symm_s([frac(14, 48), 1]), symm_s_t([-frac(1, 48), 3, 1])
            )
        else:
            assert index == 2
            self.degree = 5
            self.weights, self.points = concat(
                symm_s([frac(889, 2880), 1], [frac(5, 2880), 3]),
                symm_s_t([-frac(98, 2880), 3, 1], [frac(11, 2880), 5, 1]),
            )

        self.weights *= 4
        return
