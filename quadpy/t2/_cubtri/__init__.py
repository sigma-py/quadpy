import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["D.P. Laurie"],
    title="Algorithm 584: CUBTRI: Automatic Cubature over a Triangle",
    journal="ACM Trans. Math. Softw.",
    month="jun",
    year="1982",
    url="https://doi.org/10.1145/355993.356001",
)

this_dir = pathlib.Path(__file__).resolve().parent


def cubtri():
    """
    Se also
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    """
    return _read(this_dir / "cubtri.json", source)


register([cubtri])
