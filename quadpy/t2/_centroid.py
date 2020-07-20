from ._helpers import T2Scheme, s3


def centroid():
    weights, points = s3(1)
    return T2Scheme("Centroid rule", weights, points, 1, tol=7.850e-17)
