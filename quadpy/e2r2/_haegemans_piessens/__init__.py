import pathlib

from ...helpers import article
from .._helpers import _read, register

_source = article(
    authors=["Ann Haegemans", "Robert Piessens"],
    title="Construction of Cubature Formulas of Degree Seven and Nine Symmetric Planar Regions Using Orthogonal Polynomials",
    journal="SIAM Journal on Numerical Analysis",
    volume="14",
    number="3",
    month="jun",
    year="1977",
    pages="492-508",
    url="https://www.jstor.org/stable/2156699",
)

this_dir = pathlib.Path(__file__).resolve().parent


def haegemans_piessens_a():
    return _read(this_dir / "haegemans_piessens_a.json", _source)


def haegemans_piessens_b():
    return _read(this_dir / "haegemans_piessens_b.json", _source)


register([haegemans_piessens_a, haegemans_piessens_b])
