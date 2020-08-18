import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["Linbo Zhang", "Tao Cui", "Hui Liu"],
    title="A set of symmetric quadrature rules on triangles and tetrahedra",
    journal="Journal of Computational Mathematics",
    volume="27",
    number="1",
    month="jan",
    year="2009",
    pages="89-96",
    url="https://www.jstor.org/stable/43693493",
)

this_dir = pathlib.Path(__file__).resolve().parent


def zhang_cui_liu_1():
    return _read(this_dir / "zhang_cui_liu_1.json", source)


def zhang_cui_liu_2():
    return _read(this_dir / "zhang_cui_liu_2.json", source)


register([zhang_cui_liu_1, zhang_cui_liu_2])
