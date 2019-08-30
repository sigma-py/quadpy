import math

import numpy

from ..helpers import plot_disks


class E2rScheme:
    def __init__(self, name, weights, points, degree, citation):
        self.name = name
        self.citation = citation
        self.degree = degree

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int64]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int64]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points
        return

    def show(self, *args, **kwargs):
        import matplotlib.pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()
        return

    def plot(self, show_axes=True):
        import matplotlib.pyplot as plt

        ax = plt.gca()
        plt.axis("equal")

        if not show_axes:
            ax.set_axis_off()

        n = 2
        I0 = 2 * math.factorial(n - 1) * math.pi ** (0.5 * n) / math.gamma(0.5 * n)

        plot_disks(plt, self.points, self.weights, I0)
        return

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        return dot(f(flt(self.points).T), flt(self.weights))

    def integrate_discrete(self, data, dot=numpy.dot):
        """Quadrature where `data` are pointwise values expected to be
        defined at self.points.
        """
        flt = numpy.vectorize(float)
        return dot(data.T, flt(self.weights))

def _s8(a, b):
    return numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )


def _s4(a):
    return numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])


def _s40(a):
    return numpy.array([[+a, 0.0], [-a, 0.0], [0.0, +a], [0.0, -a]])
