import numpy as np

from ._helpers import C2Scheme


def product(scheme1d):
    schemes = scheme1d if isinstance(scheme1d, list) else 2 * [scheme1d]

    weights = np.outer(schemes[0].weights, schemes[1].weights).flatten()
    assert len(weights) > 0
    points = np.array(np.meshgrid(schemes[0].points, schemes[1].points)).reshape(2, -1)
    weights /= 4
    degree = min(s.degree for s in schemes)
    return C2Scheme(
        f"Product scheme ({scheme1d.name})",
        {"plain": np.vstack([weights, points])},
        degree,
    )
