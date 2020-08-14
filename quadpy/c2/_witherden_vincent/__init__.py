import pathlib

from ...helpers import article
from .._helpers import _read, register

_source = article(
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
    return _read(this_dir / "witherden_vincent_01.json", _source)


def witherden_vincent_03():
    return _read(this_dir / "witherden_vincent_03.json", _source)


def witherden_vincent_05():
    return _read(this_dir / "witherden_vincent_05.json", _source)


def witherden_vincent_07():
    return _read(this_dir / "witherden_vincent_07.json", _source)


def witherden_vincent_09():
    return _read(this_dir / "witherden_vincent_09.json", _source)


def witherden_vincent_11():
    return _read(this_dir / "witherden_vincent_11.json", _source)


def witherden_vincent_13():
    return _read(this_dir / "witherden_vincent_13.json", _source)


def witherden_vincent_15():
    return _read(this_dir / "witherden_vincent_15.json", _source)


def witherden_vincent_17():
    return _read(this_dir / "witherden_vincent_17.json", _source)


def witherden_vincent_19():
    return _read(this_dir / "witherden_vincent_19.json", _source)


def witherden_vincent_21():
    return _read(this_dir / "witherden_vincent_21.json", _source)


register(
    [
        witherden_vincent_01,
        witherden_vincent_03,
        witherden_vincent_05,
        witherden_vincent_07,
        witherden_vincent_09,
        witherden_vincent_11,
        witherden_vincent_13,
        witherden_vincent_15,
        witherden_vincent_17,
        witherden_vincent_19,
        witherden_vincent_21,
    ]
)
