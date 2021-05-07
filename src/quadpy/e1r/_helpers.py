import numpy as np

from ..helpers import QuadratureScheme, plot_disks_1d


class E1rScheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree):
        self.domain = "E1r"
        super().__init__(name, weights, points, degree, source=None)

    def integrate(self, f, dot=np.dot):
        x = np.array([self.points.T])
        fx = np.asarray(f(x))
        return dot(fx, self.weights)

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        if not show_axes:
            plt.gca().set_axis_off()

        plt.axis("equal")
        m = 1.1 * np.max(self.points)
        plt.plot([0, m], [0, 0], color="k")
        pts = np.column_stack([self.points, np.zeros(len(self.points))])
        plot_disks_1d(plt, pts, self.weights, total_area=1.0)
