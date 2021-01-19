import numpy as np

from ._helpers import C1Scheme


def trapezoidal():
    weights = np.array([1.0, 1.0])
    points = np.array([-1.0, 1.0])
    return C1Scheme("Trapezoidal rule", 1, weights, points)
