import numpy

from ..helpers import QuadratureScheme, plot_disks


class S2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree: int, source=None, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "S2"

    def plot(self, show_axes=False):
        from matplotlib import pyplot as plt

        ax = plt.gca()
        # change default range so that new disks will work
        plt.axis("equal")
        ax.set_xlim((-1.5, 1.5))
        ax.set_ylim((-1.5, 1.5))

        if not show_axes:
            ax.set_axis_off()

        disk1 = plt.Circle((0, 0), 1, color="k", fill=False)
        ax.add_artist(disk1)

        plot_disks(plt, self.points, self.weights, numpy.pi)
        return

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.array(f((rr + center).T))
        return numpy.pi * numpy.array(radius) ** 2 * dot(ff, self.weights)


def _z():
    return numpy.array([[0.0, 0.0]])


def _s8(a, b):
    return numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )


def _s4(a):
    return numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])


def _s40(a):
    return numpy.array([[+a, 0.0], [-a, 0.0], [0.0, +a], [0.0, -a]])


def _pm(a, b):
    return numpy.array([[+a, +b], [-a, +b], [+a, -b], [-a, -b]])


def _pmx(x):
    return numpy.array([[+x, 0], [-x, 0]])


def _pmy(y):
    return numpy.array([[0, +y], [0, -y]])
