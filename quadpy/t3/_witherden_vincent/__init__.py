import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import T3Scheme, _read, register

source = article(
    authors=["F.D. Witherden", "P.E. Vincent"],
    title="On the identification of symmetric quadrature rules for finite element methods",
    journal="Computers & Mathematics with Applications",
    volume="69",
    number="10",
    month="may",
    year="2015",
    pages="1232–1241",
    url="https://doi.org/10.1016/j.camwa.2015.03.017",
)

this_dir = pathlib.Path(__file__).resolve().parent


def witherden_vincent_01():
    degree = 1
    d = {"s4": [[1]]}
    return T3Scheme("Witherden-Vincent 1", d, degree, source)


def witherden_vincent_02():
    degree = 2
    d = {"s31": [[frac(1, 4)], [frac(1, 4) - sqrt(5) / 20]]}
    return T3Scheme("Witherden-Vincent 2", d, degree, source)


def witherden_vincent_03():
    return _read(this_dir / "witherden_vincent_03.json", source)


def witherden_vincent_05():
    return _read(this_dir / "witherden_vincent_05.json", source)


def witherden_vincent_06():
    return _read(this_dir / "witherden_vincent_06.json", source)


def witherden_vincent_07():
    return _read(this_dir / "witherden_vincent_07.json", source)


def witherden_vincent_08():
    return _read(this_dir / "witherden_vincent_08.json", source)


def witherden_vincent_09():
    return _read(this_dir / "witherden_vincent_09.json", source)


def witherden_vincent_10():
    return _read(this_dir / "witherden_vincent_10.json", source)


register(
    [
        witherden_vincent_01,
        witherden_vincent_02,
        witherden_vincent_03,
        witherden_vincent_05,
        witherden_vincent_06,
        witherden_vincent_07,
        witherden_vincent_08,
        witherden_vincent_09,
        witherden_vincent_10,
    ]
)
