import math
import pathlib

from ...helpers import article
from .._helpers import _read

_source = article(
    authors=["R. Cools", "K.J. Kim"],
    title="A survey of known and new cubature formulas for the unit disk",
    journal="Korean J. Comput. & Appl. Math.",
    volume="7",
    month="sep",
    year="2000",
    number="3",
    pages="477-485",
    # url="https://doi.org/10.1007/BF03012263"
)

this_dir = pathlib.Path(__file__).resolve().parent


def cools_kim_1():
    return _read(this_dir / "cools_kim_1.json", _source, weight_factor=1 / math.pi)


def cools_kim_2():
    return _read(this_dir / "cools_kim_2.json", _source, weight_factor=1 / math.pi)


def cools_kim_3():
    return _read(this_dir / "cools_kim_3.json", _source, weight_factor=1 / math.pi)
