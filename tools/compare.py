"""
Compare the errors of some schemes.
"""
import numpy
from matplotlib import pyplot as plt
from matplotlib import style

import quadrature

style.use("ggplot")


def f(x):
    return numpy.exp(x[0]) * numpy.exp(x[1])


schemes = (
    [quadrature.triangle.Strang(6)]
    + [quadrature.triangle.Cubtri()]
    + [quadrature.triangle.LynessJespersen(6)]
)

sample_sizes = [0.5 ** k for k in range(10)]
errors = numpy.empty((len(schemes), len(sample_sizes)))

for i, scheme in enumerate(schemes):
    for j, a in enumerate(sample_sizes):
        triangle = numpy.array([[0.0, 0.0], [a, 0.0], [0.0, a]])
        exact_value = 1.0 + numpy.exp(a) * (a - 1.0)
        val = quadrature.triangle.integrate(f, triangle, scheme)
        errors[i][j] = abs(exact_value - val)

for scheme, err in zip(schemes, errors):
    plt.loglog(sample_sizes, err, "o-", label=scheme.name)
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))

plt.show()
