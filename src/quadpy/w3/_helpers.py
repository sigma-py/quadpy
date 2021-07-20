import numpy as np

from ..helpers import QuadratureScheme, backend_to_function


class W3Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        self.domain = "W3"
        super().__init__(name, weights, points, degree, source, tol)

    def integrate(self, f, wedge, dot=np.dot):
        wedge = np.asarray(wedge)
        flt = np.vectorize(float)
        x = _transform(flt(self.points).T, wedge)
        det = _get_detJ(flt(self.points).T, wedge)
        return dot(f(x) * abs(det), flt(self.weights))

    def show(
        self,
        wedge=np.array(
            [[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[0, 0, 1], [1, 0, 1], [0, 1, 1]]],
            dtype=float,
        ),
        backend="vtk",
    ):
        edges = np.array(
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
        edges = np.moveaxis(edges, 1, 2)

        return backend_to_function[backend](
            _transform(self.points.T, wedge).T,
            self.weights,
            self.integrate(lambda _: np.ones(1), wedge),
            edges,
        )


def _transform(xi, wedge):
    mo = np.multiply.outer
    return (
        +mo(0.5 * (1.0 - xi[0] - xi[1]) * (1.0 - xi[2]), wedge[0, 0])
        + mo(0.5 * xi[0] * (1.0 - xi[2]), wedge[0, 1])
        + mo(0.5 * xi[1] * (1.0 - xi[2]), wedge[0, 2])
        + mo(0.5 * (1.0 - xi[0] - xi[1]) * (1.0 + xi[2]), wedge[1, 0])
        + mo(0.5 * xi[0] * (1.0 + xi[2]), wedge[1, 1])
        + mo(0.5 * xi[1] * (1.0 + xi[2]), wedge[1, 2])
    ).T


def _get_detJ(xi, wedge):
    mo = np.multiply.outer
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
