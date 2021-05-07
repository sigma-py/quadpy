from ._helpers import S2Scheme, register


def midpoint():
    d = {"zero2": [[1]]}
    return S2Scheme("Midpoint", d, 1)


register([midpoint])
