import math
import pathlib

from ...helpers import article
from .._helpers import _read

_source = article(
    authors=["Robert Piessens", "Ann Haegemans"],
    title="Cubature Formulas of Degree Nine for Symmetric Planar Regions",
    journal="Mathematics of Computation",
    volume="29",
    number="11",
    month="jul",
    year="1975",
    pages="810-815",
    url="https://doi.org/10.2307/2005291",
)

this_dir = pathlib.Path(__file__).resolve().parent


def piessens_haegemans():
    return _read(this_dir / "piessens_haegemans.json", _source, weight_factor=1 / math.pi)
