import numpy

from ..helpers import backend_to_function


class W3Scheme:
    def __init__(self, name, weights, points, degree, source):
        self.name = name
        self.degree = degree
        self.source = source

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points
        return

    def integrate(self, f, wedge, dot=numpy.dot):
        flt = numpy.vectorize(float)
        x = _transform(flt(self.points).T, wedge)
        det = _get_detJ(flt(self.points).T, wedge)
        return dot(f(x) * abs(det), flt(self.weights))

    def show(
        self,
        wedge=numpy.array(
            [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
            dtype=float,
        ),
        backend="vtk",
    ):
        edges = numpy.array(
            [
                [wedge[0, 0], wedge[0, 1]],
                [wedge[0, 1], wedge[0, 2]],
                [wedge[0, 2], wedge[0, 0]],
                #
                [wedge[1, 0], wedge[1, 1]],
                [wedge[1, 1], wedge[1, 2]],
                [wedge[1, 2], wedge[1, 0]],
                #
                [wedge[0, 0], wedge[1, 0]],
                [wedge[0, 1], wedge[1, 1]],
                [wedge[0, 2], wedge[1, 2]],
            ]
        )
        edges = numpy.moveaxis(edges, 1, 2)

        backend_to_function[backend](
            _transform(self.points.T, wedge).T,
            self.weights,
            self.integrate(lambda x: numpy.ones(1), wedge),
            edges,
        )
        return


def _transform(xi, wedge):
    mo = numpy.multiply.outer
    return (
        +mo(0.5 * (1.0 - xi[0] - xi[1]) * (1.0 - xi[2]), wedge[0, 0])
        + mo(0.5 * xi[0] * (1.0 - xi[2]), wedge[0, 1])
        + mo(0.5 * xi[1] * (1.0 - xi[2]), wedge[0, 2])
        + mo(0.5 * (1.0 - xi[0] - xi[1]) * (1.0 + xi[2]), wedge[1, 0])
        + mo(0.5 * xi[0] * (1.0 + xi[2]), wedge[1, 1])
        + mo(0.5 * xi[1] * (1.0 + xi[2]), wedge[1, 2])
    ).T


def _get_detJ(xi, wedge):
    mo = numpy.multiply.outer
    J0 = (
        -mo(0.5 * (1.0 - xi[2]), wedge[0, 0])
        + mo(0.5 * (1.0 - xi[2]), wedge[0, 1])
        - mo(0.5 * (1.0 + xi[2]), wedge[1, 0])
        + mo(0.5 * (1.0 + xi[2]), wedge[1, 1])
    ).T
    J1 = (
        -mo(0.5 * (1.0 - xi[2]), wedge[0, 0])
        + mo(0.5 * (1.0 - xi[2]), wedge[0, 2])
        - mo(0.5 * (1.0 + xi[2]), wedge[1, 0])
        + mo(0.5 * (1.0 + xi[2]), wedge[1, 2])
    ).T
    J2 = (
        -mo(0.5 * (1.0 - xi[0] - xi[1]), wedge[0, 0])
        - mo(0.5 * xi[0], wedge[0, 1])
        - mo(0.5 * xi[1], wedge[0, 2])
        + mo(0.5 * (1.0 - xi[0] - xi[1]), wedge[1, 0])
        + mo(0.5 * xi[0], wedge[1, 1])
        + mo(0.5 * xi[1], wedge[1, 2])
    ).T
    det = (
        +J0[0] * J1[1] * J2[2]
        + J1[0] * J2[1] * J0[2]
        + J2[0] * J0[1] * J1[2]
        - J0[2] * J1[1] * J2[0]
        - J1[2] * J2[1] * J0[0]
        - J2[2] * J0[1] * J1[0]
    )
    return det
