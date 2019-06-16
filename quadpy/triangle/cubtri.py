# -*- coding: utf-8 -*-
#
from .helpers import s1, s2, s3, concat, TriangleScheme


def Cubtri():
    """
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    Laurie, D. P.,
    Algorithm 584: CUBTRI: Automatic Cubature over a Triangle,
    ACM Trans. Math. Softw.,
    June 1982,
    <https://dl.acm.org/citation.cfm?id=356001>.
    """
    weights, bary = concat(
        s3(0.0378610912003147),
        s2(
            [0.0376204254131829, 0.1012865073234563],
            [0.0783573522441174, 0.4701420641051151],
            [0.1162714796569659, 0.2321023267750504],
            [0.0134442673751655, 0.0294808608844396],
        ),
        s1([0.0375097224552317, 0.7384168123405100, 0.2321023267750504]),
    )
    points = bary[:, 1:]
    return TriangleScheme("CUBTRI", 8, weights, points, bary)
