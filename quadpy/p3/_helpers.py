import numpy as np

from ..helpers import QuadratureScheme, backend_to_function


class P3Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "P3"

    def integrate(self, f, pyra, dot=np.dot):
        flt = np.vectorize(float)

        xi = flt(self.points).T
        x = _transform(xi, pyra)
        det = _get_det_J(pyra, xi)

        return dot(f(x) * abs(det.T), flt(self.weights))

    def show(
        self,
        pyra=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1.0]]),
        backend="mpl",
    ):
        edges = np.array(
            [
                [pyra[0], pyra[1]],
                [pyra[1], pyra[2]],
                [pyra[2], pyra[3]],
                [pyra[3], pyra[0]],
                #
                [pyra[0], pyra[4]],
                [pyra[1], pyra[4]],
                [pyra[2], pyra[4]],
                [pyra[3], pyra[4]],
            ]
        )
        edges = np.moveaxis(edges, 1, 2)

        # vol = integrate(lambda x: 1.0, pyra, Felippa(1))
        vol = 1.0
        backend_to_function[backend](
            _transform(self.points.T, pyra).T, self.weights, vol, edges
        )


def _transform(xi, pyra):
    mo = np.multiply.outer
    return (
        +mo(0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1 - xi[2]), pyra[0])
        + mo(0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1 - xi[2]), pyra[1])
        + mo(0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1 - xi[2]), pyra[2])
        + mo(0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1 - xi[2]), pyra[3])
        + mo(0.500 * (1.0 + xi[2]), pyra[4])
    ).T


def _get_det_J(pyra, xi):
    J0 = (
        -np.multiply.outer(0.125 * (1.0 - xi[1]) * (1 - xi[2]), pyra[0])
        + np.multiply.outer(0.125 * (1.0 - xi[1]) * (1 - xi[2]), pyra[1])
        + np.multiply.outer(0.125 * (1.0 + xi[1]) * (1 - xi[2]), pyra[2])
        - np.multiply.outer(0.125 * (1.0 + xi[1]) * (1 - xi[2]), pyra[3])
    ).T
    J1 = (
        -np.multiply.outer(0.125 * (1.0 - xi[0]) * (1 - xi[2]), pyra[0])
        - np.multiply.outer(0.125 * (1.0 + xi[0]) * (1 - xi[2]), pyra[1])
        + np.multiply.outer(0.125 * (1.0 + xi[0]) * (1 - xi[2]), pyra[2])
        + np.multiply.outer(0.125 * (1.0 - xi[0]) * (1 - xi[2]), pyra[3])
    ).T
    J2 = (
        -np.multiply.outer(0.125 * (1.0 - xi[0]) * (1.0 - xi[1]), pyra[0])
        - np.multiply.outer(0.125 * (1.0 + xi[0]) * (1.0 - xi[1]), pyra[1])
        - np.multiply.outer(0.125 * (1.0 + xi[0]) * (1.0 + xi[1]), pyra[2])
        - np.multiply.outer(0.125 * (1.0 - xi[0]) * (1.0 + xi[1]), pyra[3])
        + np.multiply.outer(0.500 * np.ones(1), pyra[4])
    ).T
    det = (
        +J0[0] * J1[1] * J2[2]
        + J1[0] * J2[1] * J0[2]
        + J2[0] * J0[1] * J1[2]
        - J0[2] * J1[1] * J2[0]
        - J1[2] * J2[1] * J0[0]
        - J2[2] * J0[1] * J1[0]
    )
    return det.T


def _s4(a, z):
    return [[+a, +a, z], [-a, +a, z], [+a, -a, z], [-a, -a, z]]


def _s4_0(a, z):
    return [[+a, 0.0, z], [-a, 0.0, z], [0.0, +a, z], [0.0, -a, z]]
