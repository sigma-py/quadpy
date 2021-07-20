import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["S. Wandzura", "H. Xiao"],
    title="Symmetric quadrature rules on a triangle",
    journal="Computers & Mathematics with Applications",
    volume="45",
    number="12",
    month="jun",
    year="2003",
    pages="1829-1840",
    url="https://doi.org/10.1016/S0898-1221%2803%2990004-6",
)
# Note that in the above article, the authors present the coordinates in the symmetric
# triangle [[-0.5, -sqrt(3)/2], [-0.5, +sqrt(3)/2], [1, 0]]. These have been transformed
# to barycentric coordinates here.

this_dir = pathlib.Path(__file__).resolve().parent


def wandzura_xiao_1():
    return _read(this_dir / "wx01.json", source)


def wandzura_xiao_2():
    return _read(this_dir / "wx02.json", source)


def wandzura_xiao_3():
    return _read(this_dir / "wx03.json", source)


def wandzura_xiao_4():
    return _read(this_dir / "wx04.json", source)


def wandzura_xiao_5():
    return _read(this_dir / "wx05.json", source)


def wandzura_xiao_6():
    return _read(this_dir / "wx06.json", source)


register(
    [
        wandzura_xiao_1,
        wandzura_xiao_2,
        wandzura_xiao_3,
        wandzura_xiao_4,
        wandzura_xiao_5,
        wandzura_xiao_6,
    ]
)
