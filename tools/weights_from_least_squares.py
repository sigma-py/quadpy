# -*- coding: utf-8 -*-
#
from __future__ import print_function

import numpy
import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex


def with_monomials(degree):
    exponents = numpy.concatenate([
        quadpy.helpers.partition(d, 2)
        for d in range(degree+1)
        ])

    exact_vals = numpy.array([
        integrate_monomial_over_unit_simplex(k) for k in exponents
        ])

    def fun(x):
        k = exponents.T
        # <https://stackoverflow.com/a/46689653/353337>
        s = x.shape[1:] + k.shape[1:]
        return (
            x.reshape(x.shape[0], -1, 1)**k.reshape(k.shape[0], 1, -1)
            ).prod(0).reshape(s)

    return fun, exact_vals


if __name__ == '__main__':
    scheme = quadpy.triangle.VioreanuRokhlin(19)

    x, res, rank, sv = quadpy.triangle.compute_weights(scheme)

    print('num unknowns: {}'.format(len(x)))
    print('rank A: {}'.format(rank))
    print('residual: {}'.format(res))
    assert abs(res[0]) < 1.0e-14
    print('singular values:')
    for val in sv:
        print('  {:.15e}'.format(val))

    factor = 4  # depends on the convention of the scheme
    print()
    print('solution:')
    for val in x:
        print('  {:.15e}'.format(factor*val))
