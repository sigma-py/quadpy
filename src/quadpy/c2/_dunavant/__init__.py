# ENH some closed forms
import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import C2Scheme, _read, register
from .._tyler import tyler_2 as dunavant_03

source = article(
    authors=["D.A. Dunavant"],
    title="Economical symmetrical quadrature rules for complete polynomials over a square domain",
    journal="Numerical Methods in Engineering",
    volume="21",
    number="10",
    month="oct",
    year="1985",
    pages="1777–1784",
    url="https://doi.org/10.1002/nme.1620211004",
)

this_dir = pathlib.Path(__file__).resolve().parent


def dunavant_00():
    d = {"zero2": [[1]]}
    return C2Scheme("Dunavant 0", d, 1, source, 1.0e-100)


def dunavant_01():
    d = {"d4_aa": [[frac(1, 4)], [sqrt(frac(1, 3))]]}
    return C2Scheme("Dunavant 1", d, 3, source, 4.441e-16)


def dunavant_02():
    d = {
        "d4_a0": [[frac(10, 49)], [sqrt(frac(7, 15))]],
        "d4_aa": [[frac(9, 196)], [sqrt(frac(7, 9))]],
    }
    return C2Scheme("Dunavant 2", d, 5, source, 1.305e-15)


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
    return _read(this_dir / "dunavant_09.json", source)


def dunavant_10():
    return _read(this_dir / "dunavant_10.json", source)


register(
    [
        dunavant_00,
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
    ]
)
