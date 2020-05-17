import math

from ..helpers import article, pm, untangle
from ._helpers import S2Scheme

_source = article(
    authors=["Zhongxuan Luo", "Zhaoliang Meng"],
    title="Cubature formulas over the n-sphere",
    journal="Journal of Computational and Applied Mathematics",
    year="2007",
    volume="202",
    pages="511-522",
    url="https://doi.org/10.1016/j.cam.2006.03.004",
)


def luo_meng_1():
    data = [
        (0.26179938779915, pm([0.22985042169050, 0.39811260850906])),
        (0.26179938779915, pm([0.45970084338098, 0])),
        (0.19501407677793, pm([0.33703975180642, 0.82163211980611])),
        (0.19768500492079, pm([0.81932059025905, 0.34262064294548])),
    ]
    points, weights = untangle(data)

    weights /= math.pi
    return S2Scheme("Luo-Meng 1", weights, points, 7, _source)
