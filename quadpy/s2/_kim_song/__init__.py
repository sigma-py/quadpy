import math
import pathlib

import numpy

from ...helpers import article, fsd, untangle, z
from .._helpers import S2Scheme, _read

_source = article(
    authors=["KyoungJoong Kim", "ManSuk Song"],
    title="Symmetric quadrature formulas over a unit disk",
    journal="Korean J. Comp. & Appl. Math.",
    year="1997",
    volume="4",
    pages="179-192",
    url="https://doi.org/10.1007/BF03011388",
)

this_dir = pathlib.Path(__file__).resolve().parent


def kim_song_1():
    data = [
        (1.0, z(2)),
    ]
    points, weights = untangle(data)
    points = points.astype(numpy.float)
    return S2Scheme("Kim-Song 1", weights, points, 1, _source)


def kim_song_2():
    data = [
        (0.25, fsd(2, (math.sqrt(0.5), 1))),
    ]
    points, weights = untangle(data)
    return S2Scheme("Kim-Song 2", weights, points, 3, _source)


def kim_song_3():
    data = [
        (0.25, fsd(2, (0.5, 2))),
    ]
    points, weights = untangle(data)
    return S2Scheme("Kim-Song 3", weights, points, 3, _source)


def kim_song_4():
    data = [
        (0.732786462492640, fsd(2, (0.650115167343736, 1))),
        (0.052611700904808, fsd(2, (0.888073833977115, 2))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 4", weights, points, 5, _source, 1.047e-14)


def kim_song_5():
    data = [
        (1 / 4, z(2)),
        (3 / 32, fsd(2, (0.754344479484572, 1), (0.312459714103782, 1))),
    ]
    points, weights = untangle(data)
    return S2Scheme("Kim-Song 5", weights, points, 5, _source)


def kim_song_6():
    data = [
        (2 / 27 * math.pi, fsd(2, (0.866025403784439, 1))),
        (0.387077796006226, fsd(2, (0.322914992067400, 2))),
        (0.165609800458645, fsd(2, (0.644171310389465, 2))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 6", weights, points, 7, _source)


def kim_song_7():
    data = [
        (0.071488826617391, fsd(2, (0.952458896434417, 1))),
        (0.327176874928167, fsd(2, (0.415187657878755, 1))),
        (0.005591341512851, fsd(2, (0.834794942216211, 2))),
        (0.190570560169519, fsd(2, (0.740334457173511, 1), (0.379016937530835, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 7", weights, points, 9, _source, 7.364e-14)


# TODO find issue
def kim_song_8():
    return _read(this_dir / "kim_song_08.json", _source)


def kim_song_9():
    data = [
        (0.051310052712355, fsd(2, (0.680167267076408, 2))),
        (0.208368275231940, fsd(2, (0.232463234651158, 2))),
        (0.113628206510048, fsd(2, (0.547722557505169, 2))),
        (0.126977836503225, fsd(2, (0.652159581445885, 1), (0.174745633184644, 1))),
        (0.079067977968328, fsd(2, (0.904823085572323, 1), (0.242446615072141, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 9", weights, points, 11, _source, 1.742e-14)


def kim_song_10():
    data = [
        (0.156239869837333, fsd(2, (0.283402832348825, 1))),
        (0.069447409266977, fsd(2, (0.920575474036251, 1))),
        (0.138145910534871, fsd(2, (0.649007230578002, 1))),
        (0.043376583887155, fsd(2, (0.677355106028069, 2))),
        (0.150228334140010, fsd(2, (0.412672216763289, 2))),
        (0.018199321926118, fsd(2, (0.938694583835123, 1), (0.335767600829539, 1))),
        (0.095780705939433, fsd(2, (0.752042776803954, 1), (0.379717011170077, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 10", weights, points, 13, _source, 3.073e-13)


def kim_song_11():
    data = [
        (0.125290208564338, fsd(2, (0.252863797091293, 1))),
        (0.016712625496982, fsd(2, (0.989746802511614, 1))),
        (0.109500391126365, fsd(2, (0.577728928444958, 1))),
        (0.066237455796397, fsd(2, (0.873836956645035, 1))),
        (0.026102860184358, fsd(2, (0.689299380791136, 2))),
        (0.066000934661100, fsd(2, (0.597614304667208, 2))),
        (0.127428372681720, fsd(2, (0.375416824626170, 2))),
        (0.042523065826681, fsd(2, (0.883097111318591, 1), (0.365790800400663, 1))),
        (0.081539591616413, fsd(2, (0.707438744960070, 1), (0.293030722710664, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 11", weights, points, 15, _source, 1.910e-12)


def kim_song_12():
    return _read(this_dir / "kim_song_12.json", _source)


def kim_song_13():
    return _read(this_dir / "kim_song_13.json", _source)


def kim_song_14():
    return _read(this_dir / "kim_song_14.json", _source)


def kim_song_15():
    # ENH only the first few digits are correct in the article; quadpy-optimize improved
    #     things
    return _read(this_dir / "kim_song_15.json", _source)
