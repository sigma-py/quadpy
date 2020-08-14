import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["Gatermann, Karin"],
    title="The construction of symmetric cubature formulas for the square and the triangle",
    journal="Computing",
    month="sep",
    year="1988",
    volume="40",
    number="3",
    pages="229–240",
    url="https://doi.org/10.1007/BF02251251",
)

this_dir = pathlib.Path(__file__).resolve().parent


def gatermann():
    return _read(this_dir / "gatermann.json", source)


register([gatermann])
