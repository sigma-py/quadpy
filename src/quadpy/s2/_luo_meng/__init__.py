import pathlib

from ...helpers import article
from .._helpers import _read, register

_source = article(
    authors=["Zhongxuan Luo", "Zhaoliang Meng"],
    title="Cubature formulas over the n-sphere",
    journal="Journal of Computational and Applied Mathematics",
    year="2007",
    volume="202",
    pages="511-522",
    url="https://doi.org/10.1016/j.cam.2006.03.004",
)

this_dir = pathlib.Path(__file__).resolve().parent


def luo_meng_1():
    return _read(this_dir / "luo_meng_1.json", _source)


def luo_meng_2():
    return _read(this_dir / "luo_meng_2.json", _source)


def luo_meng_3():
    return _read(this_dir / "luo_meng_3.json", _source)


def luo_meng_4():
    return _read(this_dir / "luo_meng_4.json", _source)


def luo_meng_5():
    return _read(this_dir / "luo_meng_5.json", _source)


def luo_meng_6():
    return _read(this_dir / "luo_meng_6.json", _source)


# TODO find error
# def luo_meng_7():
#     return _read(this_dir / "luo_meng_7.json", _source)


register([luo_meng_1, luo_meng_2, luo_meng_3, luo_meng_4, luo_meng_5, luo_meng_6])
