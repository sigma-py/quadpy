# -*- coding: utf-8 -*-
#
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
                (14.0/48.0, _symm_s(1.0)),
                (-1.0/48.0, _symm_s_t(3.0, 1.0))
                ]
        else:
            assert index == 2
            self.degree = 5
            data = [
                (889.0/2880.0, _symm_s(1.0)),
                (-98.0/2880.0, _symm_s_t(3.0, 1.0)),
                (5.0/2880.0, _symm_s(3.0)),
                (11.0/2880.0, _symm_s_t(5.0, 1.0)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0
        return
