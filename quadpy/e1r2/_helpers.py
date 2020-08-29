import numpy

from ..helpers import QuadratureScheme, plot_disks_1d


class E1r2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source=None):
        self.domain = "E1r2"
        super().__init__(name, weights, points, degree, source)

    def integrate(self, f, dot=numpy.dot):
        x = numpy.array([self.points.T])
        fx = numpy.asarray(f(x))
        return dot(fx, self.weights)

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        if not show_axes:
            plt.gca().set_axis_off()

        plt.axis("equal")
        m = 1.1 * numpy.max(self.points)
        plt.plot([-m, +m], [0, 0], color="k")
        pts = numpy.column_stack([self.points, numpy.zeros(len(self.points))])
        plot_disks_1d(plt, pts, self.weights, total_area=1.0)
