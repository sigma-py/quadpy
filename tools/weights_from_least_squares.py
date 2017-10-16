import numpy
import orthopy
import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
from quadpy.triangle.helpers import _s3, _s21, _s111ab
import scipy.special


# s = quadpy.triangle.Dunavant(20)
s = quadpy.triangle.Cubtri()


# WITH MONOMIALS
# ==============
# exponents = numpy.concatenate([
#     quadpy.helpers.partition(d, 2)
#     for d in range(s.degree+1)
#     ])
#
# exact_vals = numpy.array([
#     integrate_monomial_over_unit_simplex(k) for k in exponents
#     ])
#
#
# def fun(x):
#     k = exponents.T
#     # <https://stackoverflow.com/a/46689653/353337>
#     s = x.shape[1:] + k.shape[1:]
#     return (
#         x.reshape(x.shape[0], -1, 1)**k.reshape(k.shape[0], 1, -1)
#         ).prod(0).reshape(s)


# WITH ORTHOGONAL POLYNOMIALS
# ===========================
exponents = numpy.concatenate([
    quadpy.helpers.partition(d, 2)
    for d in range(s.degree+1)
    ])
# The exact integrals of the orthogonal polynomials over the triangle are 0,
# except for the one with degree 0 for which we have sqrt(2)/2.
exact_vals = numpy.zeros(len(exponents))
exact_vals[0] = numpy.sqrt(2) / 2


def fun(x):
    def f(i, j, x, y):
        alpha1, beta1 = \
            orthopy.recurrence_coefficients.jacobi(i, 0, 0, mode='numpy')
        a = orthopy.tools.evaluate_orthogonal_polynomial(
                alpha1, beta1, (x-y)/(x+y)
                )
        # a = numpy.polyval(scipy.special.jacobi(i, 0, 0), (x-y)/(x+y))

        alpha2, beta2 = \
            orthopy.recurrence_coefficients.jacobi(j, 2*i+1, 0, mode='numpy')
        b = orthopy.tools.evaluate_orthogonal_polynomial(
                alpha2, beta2, 1-2*(x+y)
                )
        # b = numpy.polyval(scipy.special.jacobi(j, 2*i+1, 0), 1-2*(x+y))

        return (
            numpy.sqrt(2*i + 1) * a * (x+y)**i
            * numpy.sqrt(2*j + 2*i + 2) * b
            )

    return numpy.array([
        f(ij[0], ij[1], x[0], x[1])
        for ij in exponents
        ]).T


a_data = []
if 's3' in s.data:
    a_data.append(fun(_s3().T[1:]))

if 's2' in s.data:
    s2_data = _s21(s.data['s2'][1])
    a_data.append(numpy.sum(fun(s2_data[1:]), axis=1))

if 's1' in s.data:
    s1_data = _s111ab(*s.data['s1'][1:])
    a_data.append(numpy.sum(fun(s1_data[1:]), axis=1))

A = numpy.concatenate(a_data).T

print(A.shape)

x, res, rank, sv = numpy.linalg.lstsq(A, exact_vals)

print
print(x*2)
print(res)
print(rank)
print(sv)
print(numpy.dot(A, x) - exact_vals)

print('{:.15e}'.format(x[0]*2))
print('{:.15e}'.format(x[1]*2))
print('{:.15e}'.format(x[-1]*2))
