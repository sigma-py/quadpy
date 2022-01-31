"""
Compare the errors of some schemes.
"""
import numpy as np
import quadrature
from matplotlib import pyplot as plt
from matplotlib import style

style.use("ggplot")


def f(x):
    return np.exp(x[0]) * np.exp(x[1])


schemes = (
    [quadrature.triangle.Strang(6)]
    + [quadrature.triangle.Cubtri()]
    + [quadrature.triangle.LynessJespersen(6)]
)

sample_sizes = [0.5**k for k in range(10)]
errors = np.empty((len(schemes), len(sample_sizes)))

for i, scheme in enumerate(schemes):
    for j, a in enumerate(sample_sizes):
        triangle = np.array([[0.0, 0.0], [a, 0.0], [0.0, a]])
        exact_value = 1.0 + np.exp(a) * (a - 1.0)
        val = quadrature.triangle.integrate(f, triangle, scheme)
        errors[i][j] = abs(exact_value - val)

for scheme, err in zip(schemes, errors):
    plt.loglog(sample_sizes, err, "o-", label=scheme.name)
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))

plt.show()
