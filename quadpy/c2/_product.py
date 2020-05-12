import numpy

from ._helpers import C2Scheme


def product(scheme1d):
    schemes = scheme1d if isinstance(scheme1d, list) else 2 * [scheme1d]

    weights = numpy.outer(schemes[0].weights, schemes[1].weights).flatten()
    points = numpy.dstack(numpy.meshgrid(schemes[0].points, schemes[1].points)).reshape(
        -1, 2
    )
    degree = min([s.degree for s in schemes])
    return C2Scheme(f"Product scheme ({scheme1d.name})", weights, points, degree)
