import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["R. Cools", "Ann Haegemans"],
    title="Another Step Forward in Searching for Cubature Formulae with a Minimal Number of Knots for the Square",
    journal="Computing",
    volume="40",
    pages="139-146",
    year="1988",
    url="https://doi.org/10.1007/BF02247942",
)

this_dir = pathlib.Path(__file__).resolve().parent


def cools_haegemans_1988_1():
    return _read(this_dir / "cools_haegemans_1988_1.json", source)


def cools_haegemans_1988_2():
    return _read(this_dir / "cools_haegemans_1988_2.json", source)


register([cools_haegemans_1988_1, cools_haegemans_1988_2])
