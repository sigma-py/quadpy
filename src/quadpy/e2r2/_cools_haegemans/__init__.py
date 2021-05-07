import pathlib

from ...helpers import techreport
from .._helpers import _read, register

_source = techreport(
    authors=["R. Cools", "A. Haegemans"],
    title="Construction of fully symmetric cubature formulae of degree 4k-3 for fully symmetric planar regions",
    year="1985",
    institution="Dept. of Computer Science, KU Leuven",
    number="Report TW 71",
    url="https://lirias.kuleuven.be/bitstream/123456789/131870/1/TW71.pdf",
)

this_dir = pathlib.Path(__file__).resolve().parent


def cools_haegemans_9_1():
    return _read(this_dir / "cools_haegemans_09_1.json", _source)


def cools_haegemans_9_2():
    return _read(this_dir / "cools_haegemans_09_2.json", _source)


def cools_haegemans_13_1():
    return _read(this_dir / "cools_haegemans_13_1.json", _source)


register([cools_haegemans_9_1, cools_haegemans_9_2, cools_haegemans_13_1])
