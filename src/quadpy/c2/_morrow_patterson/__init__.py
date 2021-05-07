import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["C.R. Morrow", "T.N.L. Patterson"],
    title="The Construction of Algebraic Cubature Formulae by the Distribution of Nodes Along Selected Lines",
    journal="SIAM J. Numer. Anal.",
    volume="22",
    number="6",
    year="1985",
    pages="1178–1190",
    url="https://doi.org/10.1137/0722071",
)

this_dir = pathlib.Path(__file__).resolve().parent


def morrow_patterson_1():
    return _read(this_dir / "morrow_patterson_1.json", source)


def morrow_patterson_2():
    return _read(this_dir / "morrow_patterson_2.json", source)


register([morrow_patterson_1, morrow_patterson_2])
