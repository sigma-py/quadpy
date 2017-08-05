# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21

from ..helpers import untangle


class SevenPoint(object):
    def __init__(self):
        self.name = 'seven-point'
        self.degree = 3
        data = [
            (0.45, _s3()),
            (0.05, _s21(0.0)),
            (2.0/15.0, _s21(0.5)),
            ]
        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
