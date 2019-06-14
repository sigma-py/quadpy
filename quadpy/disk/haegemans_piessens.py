# -*- coding: utf-8 -*-
#
"""
Ann Haegemans and Robert Piessens,
Construction of Cubature Formulas of Degree Seven and Nine Symmetric Planar Regions
Using Orthogonal Polynomials,
SIAM Journal on Numerical Analysis, Vol. 14, No. 3 (Jun., 1977), pp. 492-508,
<https://www.jstor.org/stable/2156699>.

Abstract:
A method of constructing twelve-point cubature formulas with polynomial precision seven
and nineteen-point and eighteen-point cubature formulas with polynomial precision nine
is given for planar regions and weight functions, which are symmetric in each variable.
The nodes are computed as common zeros of a set of linearly independent orthogonal
polynomials.
"""
import numpy

from .helpers import DiskScheme


def HaegemansPiessens():
    points = numpy.array(
        [
            [+0.51403903063812387088, +0.29354430793380028577],
            [+0.51403903063812387088, -0.29354430793380028577],
            [-0.51403903063812387088, +0.29354430793380028577],
            [-0.51403903063812387088, -0.29354430793380028577],
            [+0.77204117779056310600, +0.53170030383249178943],
            [+0.77204117779056310600, -0.53170030383249178943],
            [-0.77204117779056310600, +0.53170030383249178943],
            [-0.77204117779056310600, -0.53170030383249178943],
            [+0.35302323936616416863, +0.81216903083481680870],
            [+0.35302323936616416863, -0.81216903083481680870],
            [-0.35302323936616416863, +0.81216903083481680870],
            [-0.35302323936616416863, -0.81216903083481680870],
            [+0.90412681628367164440, 0.0],
            [-0.90412681628367164440, 0.0],
            [0.0, +0.98342353303608571086],
            [0.0, -0.98342353303608571086],
            [0.0, +0.57370456740697895987],
            [0.0, -0.57370456740697895987],
            [0.0, 0.0],
        ]
    )
    weights = numpy.array(
        [
            0.26278277957889186181,
            0.26278277957889186181,
            0.26278277957889186181,
            0.26278277957889186181,
            0.92628529230400928736e-1,
            0.92628529230400928736e-1,
            0.92628529230400928736e-1,
            0.92628529230400928736e-1,
            0.13543229722409008260,
            0.13543229722409008260,
            0.13543229722409008260,
            0.13543229722409008260,
            0.13413379835073218915,
            0.13413379835073218915,
            0.34823837217490051775e-1,
            0.34823837217490051775e-1,
            0.25124094894580772473,
            0.25124094894580772473,
            0.3378210604282018145,
        ]
    )
    return DiskScheme("Haegemans-Piessens", 9, weights, points)
