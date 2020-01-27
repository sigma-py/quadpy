import numpy

from ..helpers import plot_disks_1d


class LineSegmentScheme:
    def __init__(self, name, degree, weights, points, citation=None):
        self.name = name
        self.degree = degree
        self.weights = weights
        self.points = points
        self.citation = citation

    def integrate(self, f, intervals, dim=None, dot=numpy.dot):
        iv = numpy.asarray(intervals)
        x0 = 0.5 * (1.0 - self.points)
        x1 = 0.5 * (1.0 + self.points)
        x = numpy.multiply.outer(iv[0], x0) + numpy.multiply.outer(iv[1], x1)
        fx = numpy.asarray(f(x))
        if dim is None:
            # Try to guess the dimensionality of x by comparing the shapes of x and
            # f(x).
            dim = len(x.shape)
            for a, b in zip(x.shape[::-1], fx.shape[::-1]):
                if a != b:
                    break
                dim -= 1

        # numpy.sum is slower than dot() and friends, but allows for scalar input.
        diff = iv[1] - iv[0]
        len_intervals = numpy.sqrt(
            numpy.sum(diff ** 2, axis=tuple(-d for d in range(dim)))
        )
        # The factor 0.5 is from the length of the reference line [-1, 1].
        return 0.5 * len_intervals * dot(fx, self.weights)

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
