# -*- coding: utf-8 -*-
#
import numpy


def get_symmetry_code_tri(pts):
    if len(pts) == 1:
        return '_s3()'
    elif len(pts) == 3:
        # Symmetry group [[a, a, b], [a, b, a], [b, a, a]].
        # Find the equal value `a`.
        tol = 1.0e-12
        beta = pts[0] - pts[0][0]
        ct = numpy.count_nonzero(abs(beta) < tol)
        assert ct in [1, 2], beta
        val = pts[0][0] if ct == 2 else pts[0][1]
        return '_s21({:.15e})'.format(val)

    # Symmetry group [[a, b, c], [c, a, b], ...].
    assert len(pts) == 6
    # Take the two largest value from a, b, c.
    pt0 = numpy.sort(pts[0])
    return '_s111({:.15e}, {:.15e})'.format(pt0[2], pt0[1])


def get_symmetry_code_tet(pts):
    if len(pts) == 1:
        return '_s4()'
    elif len(pts) == 4:
        # Symmetry group [[a, a, a, b], [a, a, b, a], ...].
        # Find the equal value `a`.
        tol = 1.0e-12
        beta = pts[0] - pts[0][0]
        ct = numpy.count_nonzero(abs(beta) < tol)
        if ct == 1:
            return '_s31({:.15e})'.format(pts[0][1])
        else:
            assert ct == 3
            return '_s31({:.15e})'.format(pts[0][0])
    elif len(pts) == 6:
        # Symmetry group [[a, a, b, b], [a, b, b, a], ...].
        return '_s22({:.15e})'.format(pts[0][0])
    elif len(pts) == 12:
        # Symmetry group [[a, a, b, c], [a, c, a, b], ...].
        # Find the double.
        pt0 = numpy.sort(pts[0])
        tol = 1.0e-12
        for k in range(4):
            beta = pt0 - pt0[k]
            ct = numpy.count_nonzero(abs(beta) < tol)
            if ct == 2:
                return '_s211({:.15e}, {:.15e})'.format(pt0[k], pt0[k-1])

    assert len(pts) == 24
    pt0 = numpy.sort(pts[0])
    return '_s1111({:.15e}, {:.15e}, {:.15e})'.format(
            pt0[-1], pt0[-2], pt0[-3]
            )
