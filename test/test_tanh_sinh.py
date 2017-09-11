# -*- coding: utf-8 -*-
#
import pytest
import quadpy


@pytest.mark.parametrize(
    'f, f_derivatives, a, b, exact',
    # test problem 1 from
    # <http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf>
    [(lambda t: 1, {1: lambda t: 0, 2: lambda t: 0}, -1, +1, 2)]
    + [(lambda t: 1, {1: lambda t: 0, 2: lambda t: 0}, 0, +1, 1)]
    )
def test_tanh_sinh(f, f_derivatives, a, b, exact):
    # test fine error estimate
    tol = 1.0e-50
    value, error_estimate = quadpy.line_segment.tanh_sinh_quadrature(
                f, a, b, tol,
                f_derivatives=f_derivatives
                )
    assert abs(value - exact) < tol

    # test crude error estimate
    tol = 1.0e-50
    value, error_estimate = \
        quadpy.line_segment.tanh_sinh_quadrature(f, a, b, tol)
    assert abs(value - exact) < tol
    return


if __name__ == '__main__':
    test_tanh_sinh(
        lambda t: 1, {1: lambda t: 0, 2: lambda t: 0}, 2
        )
