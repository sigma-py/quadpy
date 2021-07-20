import pathlib

from ...helpers import article
from .._helpers import _read, register
from .._mysovskih import mysovskih_2, mysovskih_3

_source = article(
    authors=["Philip Rabinowitz", "Nira Richter"],
    title="Perfectly Symmetric Two-Dimensional Integration Formulas with Minimal Numbers of Points",
    journal="Mathematics of Computation",
    volume="23",
    number="108",
    month="oct",
    year="1969",
    pages="765-779",
    url="https://doi.org/10.1090/S0025-5718-1969-0258281-4",
)

this_dir = pathlib.Path(__file__).resolve().parent


def rabinowitz_richter_1():
    return _read(this_dir / "rabinowitz_richter_1.json", _source)


def rabinowitz_richter_2():
    return _read(this_dir / "rabinowitz_richter_2.json", _source)


def rabinowitz_richter_3():
    # ENH Given analytically by Mysovskih.
    return mysovskih_2()


def rabinowitz_richter_4():
    return _read(this_dir / "rabinowitz_richter_4.json", _source)


def rabinowitz_richter_5():
    return _read(this_dir / "rabinowitz_richter_5.json", _source)


def rabinowitz_richter_6():
    # ENH Given mostly analytically by Mysovskih.
    return mysovskih_3()


register(
    [
        rabinowitz_richter_1,
        rabinowitz_richter_2,
        rabinowitz_richter_3,
        rabinowitz_richter_4,
        rabinowitz_richter_5,
        rabinowitz_richter_6,
    ]
)
