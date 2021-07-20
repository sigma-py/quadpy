import pathlib

from ...helpers import article, techreport
from .._helpers import _read, register

source = techreport(
    authors=["J. Berntsen", "T.O. Espelid"],
    title="Degree 13 symmetric quadrature rules for the triangle",
    # Reports in Informatics,
    institution="Dept. of Informatics, University of Bergen",
    year="1990",
)

# This first scheme was published separately as
c2 = article(
    authors=["J. Berntsen", "T.O. Espelid"],
    title="Algorithm 706: DCUTRI: An Algorithm for Adaptive Cubature over a Collection of Triangles",
    journal="ACM Trans. Math. Softw.",
    month="sep",
    year="1992",
    url="https://doi.org/10.1145/131766.131772",
)

this_dir = pathlib.Path(__file__).resolve().parent


def dcutri():
    out = berntsen_espelid_1()
    out.source = c2
    return out


def berntsen_espelid_1():
    return _read(this_dir / "berntsen_espelid_1.json", source)


def berntsen_espelid_2():
    return _read(this_dir / "berntsen_espelid_2.json", source)


def berntsen_espelid_3():
    return _read(this_dir / "berntsen_espelid_3.json", source)


def berntsen_espelid_4():
    return _read(this_dir / "berntsen_espelid_4.json", source)


register(
    [
        dcutri,
        berntsen_espelid_1,
        berntsen_espelid_2,
        berntsen_espelid_3,
        berntsen_espelid_4,
    ]
)
