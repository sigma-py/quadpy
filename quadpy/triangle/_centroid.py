from ._helpers import TriangleScheme, s3


def centroid():
    weights, points = s3(1)
    return TriangleScheme("Centroid rule", weights, points, 1)
