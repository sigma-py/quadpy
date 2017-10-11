import numpy
import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
from quadpy.triangle.helpers import _s3, _s21, _s111ab


s = quadpy.triangle.Dunavant(20)

exponents = numpy.concatenate([
    quadpy.helpers.partition(d, 2)
    for d in range(s.degree)
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


s21_data = _s21(s.data['s21'][1])
s111_data = _s111ab(*s.data['s111'][1:])

A = numpy.concatenate([
    fun(_s3().T[1:]),
    numpy.sum(fun(s21_data[1:]), axis=0),
    numpy.sum(fun(s111_data[1:]), axis=0),
    ]).T


x, res, rank, sv = numpy.linalg.lstsq(A, exact_vals)

print
print(x*2)
print(res)
print(rank)
print(sv)
print(numpy.dot(A, x) - exact_vals)

print('{:.15e}'.format(x[0]*2))
print('{:.15e}'.format(x[-1]*2))
