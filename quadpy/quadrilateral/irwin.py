# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import _symm_s, _symm_s_t
from ..helpers import untangle


class Irwin(object):
    '''
    Joseph Oscar Irwin,
    On quadrature and cubature,
    Cambridge University Press, 1923,
    <https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ&redir_esc=y>
    '''
    def __init__(self, index):
        self.name = 'Irwin({})'.format(index)
        if index == 1:
            self.degree = 3
            data = [
                (fr(14, 48), _symm_s(1)),
                (-fr(1, 48), _symm_s_t(3, 1))
                ]
        else:
            assert index == 2
            self.degree = 5
            data = [
                (fr(889, 2880), _symm_s(1)),
                (-fr(98, 2880), _symm_s_t(3, 1)),
                (fr(5, 2880), _symm_s(3)),
                (fr(11, 2880), _symm_s_t(5, 1)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
