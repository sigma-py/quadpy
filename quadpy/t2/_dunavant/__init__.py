import pathlib

from sympy import Rational as frac

from ...helpers import article
from .._helpers import T2Scheme, _read, register

source = article(
    authors=["D.A. Dunavant"],
    title="High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle",
    journal="Article in International Journal for Numerical Methods in Engineering",
    volume="21",
    number="6",
    pages="1129-1148",
    month="jun",
    year="1985",
    url="https://doi.org/10.1002/nme.1620210612",
)

this_dir = pathlib.Path(__file__).resolve().parent


def dunavant_01():
    d = {"centroid": [[1]]}
    return T2Scheme("Dunavant 1", d, 1, source, 7.850e-17)


def dunavant_02():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Dunavant 2", d, 2, source, 2.220e-16)


def dunavant_03():
    d = {"centroid": [[-frac(9, 16)]], "d3_aa": [[frac(25, 48)], [frac(1, 5)]]}
    return T2Scheme("Dunavant 3", d, 3, source, 6.661e-16)


def dunavant_04():
    return _read(this_dir / "dunavant_04.json", source)


def dunavant_05():
    return _read(this_dir / "dunavant_05.json", source)


def dunavant_06():
    return _read(this_dir / "dunavant_06.json", source)


def dunavant_07():
    return _read(this_dir / "dunavant_07.json", source)


def dunavant_08():
    return _read(this_dir / "dunavant_08.json", source)


def dunavant_09():
    # DUP equals TRIEX 19
    return _read(this_dir / "dunavant_09.json", source)


def dunavant_10():
    return _read(this_dir / "dunavant_10.json", source)


def dunavant_11():
    return _read(this_dir / "dunavant_11.json", source)


def dunavant_12():
    return _read(this_dir / "dunavant_12.json", source)


def dunavant_13():
    return _read(this_dir / "dunavant_13.json", source)


def dunavant_14():
    return _read(this_dir / "dunavant_14.json", source)


def dunavant_15():
    return _read(this_dir / "dunavant_15.json", source)


def dunavant_16():
    return _read(this_dir / "dunavant_16.json", source)


def dunavant_17():
    return _read(this_dir / "dunavant_17.json", source)


def dunavant_18():
    return _read(this_dir / "dunavant_18.json", source)


def dunavant_19():
    return _read(this_dir / "dunavant_19.json", source)


def dunavant_20():
    return _read(this_dir / "dunavant_20.json", source)


register(
    [
        dunavant_01,
        dunavant_02,
        dunavant_03,
        dunavant_04,
        dunavant_05,
        dunavant_06,
        dunavant_07,
        dunavant_08,
        dunavant_09,
        dunavant_10,
        dunavant_11,
        dunavant_12,
        dunavant_13,
        dunavant_14,
        dunavant_15,
        dunavant_16,
        dunavant_17,
        dunavant_18,
        dunavant_19,
        dunavant_20,
    ]
)
