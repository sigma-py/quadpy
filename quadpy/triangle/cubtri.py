# -*- coding: utf-8 -*-
#
from .helpers import untangle2


class Cubtri(object):
    """
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    Laurie, D. P.,
    Algorithm 584: CUBTRI: Automatic Cubature over a Triangle,
    ACM Trans. Math. Softw.,
    June 1982,
    <http://dl.acm.org/citation.cfm?id=356001>.
    """

    def __init__(self):
        self.name = "CUBTRI"
        self.degree = 8
        self.data = {
            "s3": [0.0378610912003147],
            "s2": [
                [0.0376204254131829, 0.1012865073234563],
                [0.0783573522441174, 0.4701420641051151],
                [0.1162714796569659, 0.2321023267750504],
                [0.0134442673751655, 0.0294808608844396],
            ],
            "s1": [[0.0375097224552317, 0.7384168123405100, 0.2321023267750504]],
        }
        self.bary, self.weights = untangle2(self.data)
        self.points = self.bary[:, 1:]
        return
