import numpy as np

from ._helpers import C3Scheme


def product(scheme1d):
    schemes = scheme1d if isinstance(scheme1d, list) else 3 * [scheme1d]

    wy, wz, wx = np.meshgrid(schemes[0].weights, schemes[1].weights, schemes[2].weights)
    weights = np.vstack([wx.flatten(), wy.flatten(), wz.flatten()]).T
    weights = np.prod(weights, axis=1)
    assert len(weights) > 0
    # the order, yeah...
    y, z, x = np.meshgrid(schemes[0].points, schemes[1].points, schemes[2].points)
    points = np.array([x.flatten(), y.flatten(), z.flatten()])

    degree = min(s.degree for s in schemes)
    weights /= 8

    d = {"plain": [weights, points[0], points[1], points[2]]}
    return C3Scheme(f"Product scheme ({scheme1d.name})", d, degree)
