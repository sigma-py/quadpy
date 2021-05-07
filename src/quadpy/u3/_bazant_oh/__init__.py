import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["P. Bažant", "B.H. Oh"],
    title="Efficient Numerical Integration on the Surface of a Sphere",
    journal="ZAMM",
    volume="66",
    number="1",
    year="1986",
    pages="37-49",
    url="https://doi.org/10.1002/zamm.19860660108",
)

this_dir = pathlib.Path(__file__).resolve().parent


def bazant_oh_09():
    return _read(this_dir / "bazant_oh_09.json", source)


def bazant_oh_11():
    return _read(this_dir / "bazant_oh_11.json", source)


def bazant_oh_13():
    return _read(this_dir / "bazant_oh_13.json", source)


register([bazant_oh_09, bazant_oh_11, bazant_oh_13])
