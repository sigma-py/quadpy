# -*- coding: utf-8 -*-
#
import mpmath
import quadpy


def test_tanh_sinh():
    # test problem 1 from
    # <http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf>
    def f(t):
        # return t**2
        return 1
        # return t * mpmath.log(1 + t)

    value, error_estimate = \
        quadpy.line_segment.tanh_sinh_quadrature(f, 0, 1, 1.0e-50)
    print(value)
    print(error_estimate)
    return


if __name__ == '__main__':
    test_tanh_sinh()
