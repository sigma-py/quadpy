import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["J.I. Maeztu", "E. Sainz de la Maza"],
    title="An invariant quadrature rule of degree 11 for the tetrahedron",
    journal="C. R. Acad. Sci. Paris",
    volume="321",
    year="1995",
    pages="1263-1267",
)

this_dir = pathlib.Path(__file__).resolve().parent


def maeztu_sainz():
    # The article claims degree 11, but tests show only degree 1. :/
    # TODO find out what's going on
    return _read(this_dir / "maeztu_sainz.json", source)


register([maeztu_sainz])
