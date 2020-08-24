from ._helpers import S3Scheme, register


def midpoint():
    d = {"zero3": [[1]]}
    return S3Scheme("Midpoint", d, 1)


register([midpoint])
