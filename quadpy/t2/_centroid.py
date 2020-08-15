from ._helpers import T2Scheme, register


def centroid():
    d = {"centroid": [[1]]}
    return T2Scheme("Centroid rule", d, 1, tol=7.850e-17)


register([centroid])
