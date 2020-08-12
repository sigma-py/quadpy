import pathlib

from ...helpers import article
from .._helpers import _read

_source = article(
    authors=["Philip Rabinowitz", "Nira Richter"],
    title="Perfectly Symmetric Two-Dimensional Integration Formulas with Minimal Numbers of Points",
    journal="Mathematics of Computation",
    volume="23",
    number="108",
    month="oct",
    year="1969",
    pages="765-779",
    url="https://doi.org/10.2307/2004962",
)

this_dir = pathlib.Path(__file__).resolve().parent


def rabinowitz_richter_1():
    return _read(this_dir / "rabinowitz_richter_1.json", _source)


def rabinowitz_richter_2():
    return _read(this_dir / "rabinowitz_richter_2.json", _source)


def rabinowitz_richter_3():
    return _read(this_dir / "rabinowitz_richter_3.json", _source)


# ERR There's a misprint here somewhere.
# The sum of the weights is 406.4177935852262, not 6.283185307179584 as it should be.
# When replacing -.1010440929995067e+1 by -.1010440929995067e+3, this is fixed, but the
# scheme still only has order 3.
# We can better tackle this once we have orthogonal polynomials for E2r.
# TODO find out what's going wrong
# def rabinowitz_richter_4():
#    return _read(this_dir / "rabinowitz_richter_4.json", _source)


def rabinowitz_richter_5():
    return _read(this_dir / "rabinowitz_richter_5.json", _source)


RabinowitzRichter = {
    1: rabinowitz_richter_1,
    2: rabinowitz_richter_2,
    3: rabinowitz_richter_3,
    # 4: rabinowitz_richter_4,
    5: rabinowitz_richter_5,
}
