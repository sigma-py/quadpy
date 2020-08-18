import pathlib

from ...helpers import techreport
from .._helpers import _read, register

source = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)

this_dir = pathlib.Path(__file__).resolve().parent


def walkington_p5():
    return _read(this_dir / "walkington_p5.json", source)


register([walkington_p5])
