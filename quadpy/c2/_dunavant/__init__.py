import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import C2Scheme, _read, expand_symmetries, register

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
    d = {"zero": [[1]]}
    points, weights = expand_symmetries(d)
    return C2Scheme("Dunavant 0", weights, points, 1, source)


def dunavant_01():
    d = {"c4_aa": [[frac(1, 4)], [sqrt(frac(1, 3))]]}
    points, weights = expand_symmetries(d)
    return C2Scheme("Dunavant 1", weights, points, 3, source)


def dunavant_02():
    d = {
        "c4_a0": [[frac(40, 49)], [sqrt(frac(7, 15))]],
        "c4_aa": [[frac(9, 49)], [sqrt(frac(7, 9))]],
    }
    points, weights = expand_symmetries(d)
    weights /= 4
    return C2Scheme("Dunavant 2", weights, points, 5, source)


def dunavant_03():
    d = {
        "c4_a0": [[frac(98, 405)], [sqrt(frac(6, 7))]],
        "c4_aa": [
            [0.237431774690630, 0.520592916667394],
            [0.805979782918599, 0.380554433208316],
        ],
    }
    points, weights = expand_symmetries(d)
    weights /= 4
    return C2Scheme("Dunavant 3", weights, points, 7, source)


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
