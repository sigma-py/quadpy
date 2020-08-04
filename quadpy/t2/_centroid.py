from ._helpers import T2Scheme, expand_symmetries


def centroid():
    d = {"s3": [[1]]}
    points, weights = expand_symmetries(d)
    return T2Scheme("Centroid rule", weights, points, 1, tol=7.850e-17)
