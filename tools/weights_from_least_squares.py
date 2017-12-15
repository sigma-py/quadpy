# -*- coding: utf-8 -*-
#
from __future__ import print_function

import numpy
import orthopy
import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex

from orthopy.triangle.helpers import _s3, _s21, _s111ab, _rot_ab


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


def compute_weights(scheme):
    '''Computes weights from points for a given scheme. Useful for
    cross-checking scheme integrity.
    '''
    def eval_orthopolys(bary):
        return numpy.concatenate(
            orthopy.triangle.orth_tree(scheme.degree, bary, 'normal')
            )

    fun = eval_orthopolys

    a_data = []
    if 's3' in scheme.data:
        a_data.append(fun(_s3().T))

    if 's2' in scheme.data:
        d = numpy.array(scheme.data['s2']).T
        s2_data = _s21(d[1])
        a_data.append(numpy.sum(fun(s2_data), axis=1))

    if 's1' in scheme.data:
        d = numpy.array(scheme.data['s1']).T
        s1_data = _s111ab(*d[1:])
        a_data.append(numpy.sum(fun(s1_data), axis=1))

    if 'rot' in scheme.data:
        d = numpy.array(scheme.data['rot']).T
        rot_data = _rot_ab(*d[1:])
        a_data.append(numpy.sum(fun(rot_data), axis=1))

    A = numpy.column_stack(a_data)

    exact_vals = numpy.zeros(len(A))
    exact_vals[0] = numpy.sqrt(2) / 2

    return numpy.linalg.lstsq(A, exact_vals)


if __name__ == '__main__':
    scheme = quadpy.triangle.Papanicolopulos('rot', 8)
    factor = 2  # depends on the convention of the scheme

    x, res, rank, sv = quadpy.triangle.compute_weights(scheme)

    print('num unknowns: {}'.format(len(x)))
    print('rank A: {}'.format(rank))
    print('residual: {}'.format(res))
    assert abs(res[0]) < 1.0e-14
    print('singular values:')
    for val in sv:
        print('  {:.15e}'.format(val))

    print()
    print('solution:')
    for val in x:
        print('  {:.15e}'.format(factor*val))
