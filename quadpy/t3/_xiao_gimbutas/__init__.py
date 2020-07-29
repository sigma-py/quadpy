# TODO all points are "plain" -- identiy symmetries and list them as such
import pathlib

from ...helpers import article
from .._helpers import _read

source = article(
    authors=["Hong Xiao", "Zydrunas Gimbutas"],
    title="A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions",
    journal="Computers & Mathematics with Applications",
    volume="59",
    number="2",
    month="jan",
    year="2010",
    pages="663–676",
    url="https://doi.org/10.1016/j.camwa.2009.10.027",
)
# Data extracted from
# https://people.sc.fsu.edu/~jburkardt/f_src/triangle_symq_rule/triangle_symq_rule.f90

this_dir = pathlib.Path(__file__).resolve().parent


def xiao_gimbutas_01():
    return _read(this_dir / "xg01.json", source)


def xiao_gimbutas_02():
    return _read(this_dir / "xg02.json", source)


def xiao_gimbutas_03():
    return _read(this_dir / "xg03.json", source)


def xiao_gimbutas_04():
    return _read(this_dir / "xg04.json", source)


def xiao_gimbutas_05():
    return _read(this_dir / "xg05.json", source)


def xiao_gimbutas_06():
    return _read(this_dir / "xg06.json", source)


def xiao_gimbutas_07():
    return _read(this_dir / "xg07.json", source)


def xiao_gimbutas_08():
    return _read(this_dir / "xg08.json", source)


def xiao_gimbutas_09():
    return _read(this_dir / "xg09.json", source)


def xiao_gimbutas_10():
    return _read(this_dir / "xg10.json", source)


def xiao_gimbutas_11():
    return _read(this_dir / "xg11.json", source)


def xiao_gimbutas_12():
    return _read(this_dir / "xg12.json", source)


def xiao_gimbutas_13():
    return _read(this_dir / "xg13.json", source)


def xiao_gimbutas_14():
    return _read(this_dir / "xg14.json", source)


def xiao_gimbutas_15():
    return _read(this_dir / "xg15.json", source)
