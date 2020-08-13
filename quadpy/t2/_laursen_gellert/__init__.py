import pathlib

import numpy
from sympy import Rational as frac

from ...helpers import article
from .._helpers import T2Scheme, _read

source = article(
    authors=["M.E. Laursen", "M. Gellert"],
    title="Some criteria for numerically integrated matrices and quadrature formulas for triangles",
    journal="International Journal for Numerical Methods in Engineering",
    volume="12",
    number="1",
    year="1978",
    pages="67–76",
    url="https://doi.org/10.1002/nme.1620120107",
)

this_dir = pathlib.Path(__file__).resolve().parent


def laursen_gellert_01():
    d = {"s3": [[1]]}
    return T2Scheme("Laursen-Gellert 1", d, 1, source)


def laursen_gellert_02a():
    d = {"s2": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Laursen-Gellert 2a", d, 2, source)


def laursen_gellert_02b():
    d = {"s2": [[frac(1, 3)], [frac(1, 2)]]}
    return T2Scheme("Laursen-Gellert 2b", d, 2, source)


def laursen_gellert_03():
    d = {"s3": [[-frac(9, 16)]], "s2": [[frac(25, 48)], [frac(1, 5)]]}
    return T2Scheme("Laursen-Gellert 3", d, 3, source)


def laursen_gellert_04():
    roots = numpy.polynomial.polynomial.polyroots([-1, 15, -60, 60])
    d = {"s1": [[1 / 6], [roots[2]], [roots[1]]]}
    return T2Scheme("Laursen-Gellert 4", d, 3, source)


def laursen_gellert_05():
    return _read(this_dir / "laursen_gellert_05.json", source)


def laursen_gellert_06():
    d = {"s3": [[3 / 8]], "s1": [[5 / 48], [0.736712498968435], [0.237932366472434]]}
    return T2Scheme("Laursen-Gellert 6", d, 4, source)


def laursen_gellert_07():
    d = {
        "s3": [[9 / 40]],
        "s2": [
            [0.125939180544827, 0.132394152788506],
            [0.101286507323456, 0.470142064105115],
        ],
    }
    return T2Scheme("Laursen-Gellert 7", d, 5, source)


def laursen_gellert_08():
    return _read(this_dir / "laursen_gellert_08.json", source)


def laursen_gellert_09():
    return _read(this_dir / "laursen_gellert_09.json", source)


def laursen_gellert_10():
    return _read(this_dir / "laursen_gellert_10.json", source)


def laursen_gellert_11():
    return _read(this_dir / "laursen_gellert_11.json", source)


def laursen_gellert_12():
    return _read(this_dir / "laursen_gellert_12.json", source)


def laursen_gellert_13():
    return _read(this_dir / "laursen_gellert_13.json", source)


def laursen_gellert_14():
    return _read(this_dir / "laursen_gellert_14.json", source)


def laursen_gellert_15a():
    return _read(this_dir / "laursen_gellert_15a.json", source)


def laursen_gellert_15b():
    return _read(this_dir / "laursen_gellert_15b.json", source)
