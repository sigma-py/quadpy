import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["Karin Gatermann"],
    title="Linear Representations of Finite Groups and The Ideal Theoretical Construction of G-Invariant Cubature Formulas",
    journal="Numerical Integration",
    pages="25-35",
    note="Part of the NATO ASI Series book series (ASIC, volume 357)",
)

this_dir = pathlib.Path(__file__).resolve().parent


def gatermann():
    return _read(this_dir / "gatermann.json", source)


register([gatermann])
