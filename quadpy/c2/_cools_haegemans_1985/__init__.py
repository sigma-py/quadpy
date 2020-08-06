# TODO There are three more schemes in the technical report
import pathlib

from ...helpers import techreport
from .._helpers import _read

source = techreport(
    authors=["R. Cools", "A. Haegemans"],
    title="Construction of fully symmetric cubature formulae of degree 4k-3 for fully symmetric planar regions",
    year="1985",
    institution="Dept. of Computer Science, KU Leuven",
    number="Report TW 71",
    url="https://lirias.kuleuven.be/bitstream/123456789/131870/1/TW71.pdf",
)

this_dir = pathlib.Path(__file__).resolve().parent


def cools_haegemans_1985_9_1():
    return _read(this_dir / "cools_haegemans_1985_09_1.json", source)


def cools_haegemans_1985_13_1():
    return _read(this_dir / "cools_haegemans_1985_13_1.json", source)


def cools_haegemans_1985_13_2():
    return _read(this_dir / "cools_haegemans_1985_13_2.json", source)


def cools_haegemans_1985_13_3():
    return _read(this_dir / "cools_haegemans_1985_13_3.json", source)
