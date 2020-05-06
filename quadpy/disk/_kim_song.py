import warnings

import numpy

from ..helpers import article, untangle, z, fsd
from ._helpers import DiskScheme

_citation = article(
    authors=["KyoungJoong Kim", "ManSuk Song"],
    title="Symmetric quadrature formulas over a unit disk",
    journal="Korean J. Comp. & Appl. Math.",
    year="1997",
    volume="4",
    pages="179-192",
    url="https://doi.org/10.1007/BF03011388",
)


def kim_song_1():
    data = [
        (3.141592653589793, z(2)),
    ]
    points, weights = untangle(data)
    points = points.astype(numpy.float)
    return DiskScheme("Kim-Song 1", weights, points, 1, _citation)


def kim_song_2():
    data = [
        (0.785398163397448, fsd(2, (0.707106781186548, 1))),
    ]
    points, weights = untangle(data)
    return DiskScheme("Kim-Song 2", weights, points, 3, _citation)


def kim_song_3():
    data = [
        (0.785398163397448, fsd(2, (0.5, 2))),
    ]
    points, weights = untangle(data)
    return DiskScheme("Kim-Song 3", weights, points, 3, _citation)


def kim_song_4():
    data = [
        (0.732786462492640, fsd(2, (0.650115167343736, 1))),
        (0.052611700904808, fsd(2, (0.888073833977115, 2))),
    ]
    points, weights = untangle(data)
    return DiskScheme("Kim-Song 4", weights, points, 5, _citation)


def kim_song_5():
    data = [
        (0.785398163397449, z(2)),
        (0.294524311274043, fsd(2, (0.754344479484572, 1), (0.312459714103782, 1))),
    ]
    points, weights = untangle(data)
    return DiskScheme("Kim-Song 5", weights, points, 5, _citation)


def kim_song_6():
    data = [
        (0.232710566932577, fsd(2, (0.866025403784439, 1))),
        (0.387077796006226, fsd(2, (0.322914992067400, 2))),
        (0.165609800458645, fsd(2, (0.644171310389465, 2))),
    ]
    points, weights = untangle(data)
    return DiskScheme("Kim-Song 6", weights, points, 7, _citation)


def kim_song_7():
    data = [
        (0.071488826617391, fsd(2, (0.952458896434417, 1))),
        (0.327176874928167, fsd(2, (0.415187657878755, 1))),
        (0.005591341512851, fsd(2, (0.834794942216211, 2))),
        (0.190570560169519, fsd(2, (0.740334457173511, 1), (0.379016937530835, 1))),
    ]
    points, weights = untangle(data)
    return DiskScheme("Kim-Song 7", weights, points, 9, _citation)


def kim_song_8():
    data = [
        (-0.103220620849181, z(2)),
        (0.000561120485457, fsd(2, (1.081713912649978, 2))),
        (0.234185472691706, fsd(2, (0.329855087221954, 2))),
        (0.123731369524117, fsd(2, (0.163309595229403, 2))),
        (0.194337435615691, fsd(2, (0.582963679648472, 2))),
        (0.129193960146342, fsd(2, (0.862551409331466, 1), (0.191706749096887, 1))),
    ]
    points, weights = untangle(data)
    # ERR article claims degree 11
    warnings.warn("Kim-Song claim degree 11, but the scheme is only degree 1.")
    return DiskScheme("Kim-Song 8", weights, points, 1, _citation)
