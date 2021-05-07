import numpy as np

from ._helpers import C1Scheme


def midpoint():
    weights = np.array([2.0])
    points = np.array([0.0])
    return C1Scheme("Midpoint rule", 1, weights, points)
