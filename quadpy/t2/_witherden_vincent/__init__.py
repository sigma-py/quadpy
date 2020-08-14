import pathlib

from ...helpers import article
from .._helpers import _read, register

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
    return _read(this_dir / "wv01.json", source)


def witherden_vincent_02():
    return _read(this_dir / "wv02.json", source)


def witherden_vincent_04():
    return _read(this_dir / "wv04.json", source)


def witherden_vincent_05():
    return _read(this_dir / "wv05.json", source)


def witherden_vincent_06():
    return _read(this_dir / "wv06.json", source)


def witherden_vincent_07():
    return _read(this_dir / "wv07.json", source)


def witherden_vincent_08():
    return _read(this_dir / "wv08.json", source)


def witherden_vincent_09():
    return _read(this_dir / "wv09.json", source)


def witherden_vincent_10():
    return _read(this_dir / "wv10.json", source)


def witherden_vincent_11():
    return _read(this_dir / "wv11.json", source)


def witherden_vincent_12():
    return _read(this_dir / "wv12.json", source)


def witherden_vincent_13():
    return _read(this_dir / "wv13.json", source)


def witherden_vincent_14():
    return _read(this_dir / "wv14.json", source)


def witherden_vincent_15():
    return _read(this_dir / "wv15.json", source)


def witherden_vincent_16():
    return _read(this_dir / "wv16.json", source)


def witherden_vincent_17():
    return _read(this_dir / "wv17.json", source)


def witherden_vincent_18():
    return _read(this_dir / "wv18.json", source)


def witherden_vincent_19():
    return _read(this_dir / "wv19.json", source)


def witherden_vincent_20():
    return _read(this_dir / "wv20.json", source)


register(
    [
        witherden_vincent_01,
        witherden_vincent_02,
        witherden_vincent_04,
        witherden_vincent_05,
        witherden_vincent_06,
        witherden_vincent_07,
        witherden_vincent_08,
        witherden_vincent_09,
        witherden_vincent_10,
        witherden_vincent_11,
        witherden_vincent_12,
        witherden_vincent_13,
        witherden_vincent_14,
        witherden_vincent_15,
        witherden_vincent_16,
        witherden_vincent_17,
        witherden_vincent_18,
        witherden_vincent_19,
        witherden_vincent_20,
    ]
)
