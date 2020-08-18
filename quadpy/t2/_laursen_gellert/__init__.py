import pathlib

from sympy import Rational as frac

from ...helpers import article
from .._classical import centroid as laursen_gellert_01
from .._helpers import T2Scheme, _read, register
from .._strang_fix_cowper import strang_fix_cowper_03 as laursen_gellert_03  # DUP
from .._strang_fix_cowper import strang_fix_cowper_04 as laursen_gellert_04
from .._strang_fix_cowper import strang_fix_cowper_05 as laursen_gellert_05

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


def laursen_gellert_02a():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Laursen-Gellert 2a", d, 2, source)


def laursen_gellert_02b():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 2)]]}
    return T2Scheme("Laursen-Gellert 2b", d, 2, source)


def laursen_gellert_06():
    d = {
        "centroid": [[3 / 8]],
        "d3_ab": [[5 / 48], [0.736712498968435], [0.237932366472434]],
    }
    return T2Scheme("Laursen-Gellert 6", d, 4, source)


def laursen_gellert_07():
    d = {
        "centroid": [[9 / 40]],
        "d3_aa": [
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


register(
    [
        laursen_gellert_01,
        laursen_gellert_02a,
        laursen_gellert_02b,
        laursen_gellert_03,
        laursen_gellert_04,
        laursen_gellert_05,
        laursen_gellert_06,
        laursen_gellert_07,
        laursen_gellert_08,
        laursen_gellert_09,
        laursen_gellert_10,
        laursen_gellert_11,
        laursen_gellert_12,
        laursen_gellert_13,
        laursen_gellert_14,
        laursen_gellert_15a,
        laursen_gellert_15b,
    ]
)
