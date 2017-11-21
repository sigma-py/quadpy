# -*- coding: utf-8 -*-
#
import numpy

import helpers


def test():
    x = numpy.array([[2, 3, 5]]).T
    vals = helpers.evaluate_all_monomials(x, max_degree=3)
    tol = 1.0e-14
    assert numpy.all(abs(vals[0] - [1]) < tol)
    assert numpy.all(abs(vals[1] - [2, 3, 5]) < tol)
    assert numpy.all(abs(vals[2] - [4, 6, 10, 9, 15, 25]) < tol)
    assert numpy.all(
        abs(vals[3] - [8, 12, 20, 18, 30, 50, 27, 45, 75, 125]) < tol
        )
    return
