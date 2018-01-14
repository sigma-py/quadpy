# -*- coding: utf-8 -*-
#
from mpmath import mp

import quadpy


def test_hermite():
    scheme = quadpy.e1r2.Hermite(4, mode='mpmath', decimal_places=51)

    tol = 1.0e-50

    x1 = mp.sqrt((3 - mp.sqrt(6)) / 2)
    x2 = mp.sqrt((3 + mp.sqrt(6)) / 2)
    assert (abs(scheme.points - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = mp.sqrt(mp.pi) / 4 / (3 - mp.sqrt(6))
    w2 = mp.sqrt(mp.pi) / 4 / (3 + mp.sqrt(6))

    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()
    return


if __name__ == '__main__':
    test_hermite()
