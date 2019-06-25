# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import plot_disks


class DiskScheme:
    def __init__(self, name, weights, points, degree: int, citation=None):
        self.name = name
        self.weights = weights
        self.points = points
        self.degree = degree
        self.citation = citation
        return

    def show(self, *args, **kwargs):
        import matplotlib.pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()
        return

    def plot(self, show_axes=False):
        import matplotlib.pyplot as plt

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
        return numpy.array(radius) ** 2 * dot(ff, self.weights)


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
