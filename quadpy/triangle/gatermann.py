# -*- coding: utf-8 -*-
#
"""
Karin Gatermann,
The construction of symmetric cubature formulas for the square and the triangle,
Computing, September 1988, Volume 40, Issue 3, pp 229–240,
<https://doi.org/10.1007/BF02251251>.

Abstract:
The weights and nodes of a symmetric cubature formula are determined by solving a system
of nonlinear equations. The number of equations and their structure are investigated for
symmetric cubature formulas for the square and the triangle. A new cubature formula of
degree 7 with 12 nodes is given for the triangle.
"""
from .helpers import TriangleScheme, rot_ab


def Gatermann():
    weights, bary = rot_ab(
        [0.2651702815743450e-01, 0.6238226509439084e-01, 0.6751786707392436e-01],
        [0.4388140871444811e-01, 0.5522545665692000e-01, 0.3215024938520156],
        [0.2877504278497528e-01, 0.3432430294509488e-01, 0.6609491961867980],
        [0.6749318700980879e-01, 0.5158423343536001, 0.2777161669764050],
    )
    weights *= 2
    points = bary[:, 1:]
    return TriangleScheme("Gatermann", 7, weights, points, bary)
