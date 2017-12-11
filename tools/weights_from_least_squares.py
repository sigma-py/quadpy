# -*- coding: utf-8 -*-
#
from __future__ import print_function

import numpy
import orthopy
import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
from quadpy.triangle.helpers import _s3, _s21, _s111ab


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
    scheme = quadpy.triangle.Cubtri()

    def eval_orthopolys(x):
        bary = numpy.array([x[0], x[1], 1.0-x[0]-x[1]])
        out = numpy.concatenate(
            orthopy.triangle.orth_tree(scheme.degree+1, bary, 'normal')
            )
        return out

    fun = eval_orthopolys

    a_data = []
    if 's3' in scheme.data:
        a_data.append(fun(_s3().T[1:]))

    if 's2' in scheme.data:
        s2_data = _s21(scheme.data['s2'][1])
        a_data.append(numpy.sum(fun(s2_data[1:]), axis=1))

    if 's1' in scheme.data:
        s1_data = _s111ab(*scheme.data['s1'][1:])
        a_data.append(numpy.sum(fun(s1_data[1:]), axis=1))

    A = numpy.column_stack(a_data)

    exact_vals = numpy.zeros(len(A))
    exact_vals[0] = 1.0

    x, res, rank, sv = numpy.linalg.lstsq(A, exact_vals)

    print('A.shape = {} x {}'.format(*A.shape))
    print('rank: {}'.format(rank))
    print('residual: {}'.format(res))
    print('singular values:')
    for val in sv:
        print('  {:.15e}'.format(val))

    print()
    print('solution:')
    for val in x:
        print('  {:.15e}'.format(2*val))
