import json
import os

from ...helpers import article
from .._helpers import TriangleScheme, untangle2

citation = article(
    authors=["S. Wandzura", "H. Xiao"],
    title="Symmetric quadrature rules on a triangle",
    journal="Computers & Mathematics with Applications",
    volume="45",
    number="12",
    month="jun",
    year="2003",
    pages="1829-1840",
    url="https://doi.org/10.1016/S0898-1221(03)90004-6",
)
# Note that in the above article, the authors present the coordinates in the symmetric
# triangle [[-0.5, -sqrt(3)/2], [-0.5, +sqrt(3)/2], [1, 0]]. These have been transformed
# to barycentric coordinates here.


def _read(index):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"wx{index:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)
    degree = data.pop("degree")

    data = data
    points, weights = untangle2(data)
    return TriangleScheme(f"Wandzura-Xiao {index}", weights, points, degree, citation)


def wandzura_xiao_1():
    return _read(1)


def wandzura_xiao_2():
    return _read(2)


def wandzura_xiao_3():
    return _read(3)


def wandzura_xiao_4():
    return _read(4)


def wandzura_xiao_5():
    return _read(5)


def wandzura_xiao_6():
    return _read(6)
