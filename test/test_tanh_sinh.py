# -*- coding: utf-8 -*-
#
import quadpy


def test_tanh_sinh():
    # test problem 1 from
    # <http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf>
    def f(t):
        # return t**2
        return 1
        # return t * mpmath.log(1 + t)

    tol = 1.0e-50
    value, error_estimate = \
        quadpy.line_segment.tanh_sinh_quadrature(
                f, -1, +1, tol,
                f_derivatives={
                    1: lambda t: 0,
                    2: lambda t: 0,
                    }
                )
    assert abs(value - 2) < tol
    return


if __name__ == '__main__':
    test_tanh_sinh()
