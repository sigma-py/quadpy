# -*- coding: utf-8 -*-
#
from __future__ import division
import math

import numpy

from ..helpers import z, untangle


class Radon(object):
    '''
    J. Radon,
    Zur mechanischen Kubatur,
    Monatshefte für Mathematik (1948),
    Volume: 52, page 286-300, ISSN: 0026-9255; 1436-5081/e.
    <https://eudml.org/doc/176796>.
    '''
    def __init__(self, alpha):
        self.degree = 5

        r = math.sqrt((alpha+4) / (alpha+6))
        s = math.sqrt((alpha+4) / 4 / (alpha+6))
        t = math.sqrt(3 * (alpha+4) / 4 / (alpha+6))

        A = 4 / (alpha+4)**2
        B = (alpha+2) * (alpha+6) / 6 / (alpha+4)**2

        data = [
            (A, z(2)),
            (B, numpy.array([
                [+r, 0.0],
                [-r, 0.0],
                ])),
            # Stroud is missing +- in front of t.
            (B, numpy.array([
                [+s, +t],
                [-s, +t],
                [+s, -t],
                [-s, -t],
                ])),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= math.pi
        return
