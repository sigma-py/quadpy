# -*- coding: utf-8 -*-
#
from mpmath import mp
import pytest
import quadpy
import sympy


@pytest.mark.parametrize(
    'f, f_derivatives, a, b, exact',
    [({
        0: lambda t: 1,
        1: lambda t: 0,
        2: lambda t: 0
        }, -1, +1, 2)]
    + [({
        0: lambda t: 1,
        1: lambda t: 0,
        2: lambda t: 0
        }, 0, +1, 1)]
    + [({
        0: lambda t: t,
        1: lambda t: 1,
        2: lambda t: 0,
        }, -0, +1, sympy.Rational(1, 2))]
    + [({
        0: lambda t: t**2,
        1: lambda t: 2*t,
        2: lambda t: 2,
        }, -1, +1, mp.mpf(2)/3)]
    # Some test problems from
    # <http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf>
    # Test problem 1:
    + [({
        0: lambda t: t * mp.log(1+t),
        1: lambda t: t/(t+1) + mp.log(t+1),
        2: lambda t: (t+2) / (t+1)**2,
        }, 0, +1, 0.25)]
    )
def test_tanh_sinh(f, a, b, exact):
    # test fine error estimate
    tol = 1.0e-50
    value, error_estimate = quadpy.line_segment.tanh_sinh_quadrature(
                f[0], a, b, tol,
                f_derivatives={1: f[1], 2: f[2]}
                )
    assert abs(value - exact) < tol

    # test crude error estimate
    tol = 1.0e-50
    value, error_estimate = \
        quadpy.line_segment.tanh_sinh_quadrature(f[0], a, b, tol)
    assert abs(value - exact) < tol
    return


if __name__ == '__main__':
    # test_tanh_sinh(
    #     {
    #         0: lambda t: 1,
    #         1: lambda t: 0,
    #         2: lambda t: 0
    #     }, -0, +1, 1
    #     )
    # test_tanh_sinh(
    #     {
    #         0: lambda t: (t+1)/2,
    #         1: lambda t: sympy.Rational(1)/2,
    #         2: lambda t: 0,
    #     },
    #     -1, +1, 1
    #     )
    test_tanh_sinh(
        {
            0: lambda t: t,
            1: lambda t: 1,
            2: lambda t: 0,
        },
        0, +1, sympy.Rational(1, 2)
        )
    # test_tanh_sinh(
    #     {
    #         0: lambda t: t**2,
    #         1: lambda t: 2*t,
    #         2: lambda t: 2,
    #     },
    #     -1, +1, sympy.Rational(2, 3)
    #     )
    # test_tanh_sinh(
    #     lambda t: t * mp.log(1+t), {
    #         1: lambda t: t/(t+1) + mp.log(t+1),
    #         2: lambda t: (t+2) / (t+1)**2,
    #     }, 0, +1, 0.25
    #     )
