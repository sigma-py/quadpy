from ._helpers import T2Scheme


def centroid():
    d = {"s3": [[1]]}
    return T2Scheme("Centroid rule", d, 1, tol=7.850e-17)
