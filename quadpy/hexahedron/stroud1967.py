# -*- coding: utf-8 -*-
#
from .helpers import rss_pm, z

from ..helpers import untangle


class Stroud1967(object):
    '''
    A.H. Stroud,
    Some fifth degree integration formulas for symmetric regions II,
    Numerische Mathematik, Volume 9 Issue 5, April 1967, Pages 460-468
    <https://dx.doi.org/10.1007/BF02162160>.
    '''
    def __init__(self):
        self.degree = 5
        data = [
            (4.0/19.0, z()),
            (0.681234189096984e-1, rss_pm(0.880304406699309, -0.495848171425710)),
            (0.634555284587268e-1, rss_pm(0.252937117428425e-1, 0.795621422164093)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return
