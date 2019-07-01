import numpy

from ..helpers import plot_disks_1d


class LineSegmentScheme:
    def __init__(self, name, degree, weights, points, citation=None):
        self.name = name
        self.degree = degree
        self.weights = weights
        self.points = points
        self.citation = citation
        return

    def integrate(self, f, interval, dot=numpy.dot):
        xi = self.points
        x = +numpy.multiply.outer(0.5 * (1.0 - xi), interval[0]) + numpy.multiply.outer(
            0.5 * (1.0 + xi), interval[1]
        )
        x = x.T
        diff = interval[1] - interval[0]
        # numpy.sum is slower than dot() and friends, but allows for scalar input.
        len_intervals = numpy.sqrt(numpy.sum(diff ** 2, axis=-1))
        # The factor 0.5 is from the length of the reference line [-1, 1].
        return 0.5 * len_intervals * dot(f(x), self.weights)

    def integrate_split(self, f, a, b, n, dot=numpy.dot):
        """Integrates f between a and b with n subintervals.
        """
        # prepare the intervals
        x = numpy.linspace(a, b, n + 1)
        intervals = numpy.expand_dims(numpy.stack([x[:-1], x[1:]]), axis=-1)
        # integrate
        out = self.integrate(f, intervals, dot=dot)[0]
        # sum over the intervals
        return numpy.sum(out)

    def show(self, *args, **kwargs):
        import matplotlib.pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()
        return

    def plot(self, interval=numpy.array([[-1.0], [1.0]]), show_axes=False):
        import matplotlib.pyplot as plt

        # change default range so that new disks will work
        plt.axis("equal")
        # ax.set_xlim((-1.5, 1.5))
        # ax.set_ylim((-1.5, 1.5))

        if not show_axes:
            plt.gca().set_axis_off()

        plt.plot(interval, [0, 0], color="k")

        pts = numpy.column_stack([self.points, numpy.zeros(len(self.points))])

        total_area = interval[1] - interval[0]
        plot_disks_1d(plt, pts, self.weights, total_area)
        return
