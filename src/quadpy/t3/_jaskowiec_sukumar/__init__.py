import json
import pathlib

import numpy as np

from ...helpers import article
from .._helpers import T3Scheme, register

source = article(
    authors=["Jan Jaśkowiec", "N. Sukumar"],
    title="High-order cubature rules for tetrahedra",
    journal="Numerical Methods in Engineering",
    volume="121",
    number="11",
    year="2020",
    pages="2418-2436",
    url="https://doi.org/10.1002/nme.6313",
)
# Data from
# https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6313

this_dir = pathlib.Path(__file__).resolve().parent


def _read(string, tol=1.0e-14):
    filename = f"js{string}.json"
    with open(this_dir / filename) as f:
        data = json.load(f)

    degree = data["degree"]

    points = np.array(data["points"]).T
    points = np.array([points[0], points[1], points[2], 1.0 - np.sum(points, axis=0)])
    weights = np.array(data["weights"])

    d = {"plain": [weights, points[0], points[1], points[2], points[3]]}

    return T3Scheme(f"Jaśkowiec-Sukumar {string}", d, degree, source, tol)


def jaskowiec_sukumar_02():
    return _read("02")


def jaskowiec_sukumar_03():
    return _read("03")


def jaskowiec_sukumar_04():
    return _read("04")


def jaskowiec_sukumar_05():
    return _read("05")


def jaskowiec_sukumar_06():
    return _read("06")


def jaskowiec_sukumar_07():
    return _read("07")


def jaskowiec_sukumar_08():
    return _read("08")


def jaskowiec_sukumar_09():
    return _read("09")


def jaskowiec_sukumar_10():
    return _read("10")


def jaskowiec_sukumar_11():
    return _read("11")


def jaskowiec_sukumar_12():
    return _read("12")


def jaskowiec_sukumar_13():
    return _read("13")


def jaskowiec_sukumar_14():
    return _read("14")


def jaskowiec_sukumar_15():
    return _read("15", 3.687e-11)


def jaskowiec_sukumar_16():
    return _read("16", 2.653e-11)


def jaskowiec_sukumar_17():
    return _read("17", 3.831e-11)


def jaskowiec_sukumar_18():
    return _read("18", 3.470e-11)


def jaskowiec_sukumar_19a():
    return _read("19a", 4.679e-11)


def jaskowiec_sukumar_19b():
    return _read("19b", 4.580e-11)


def jaskowiec_sukumar_20():
    return _read("20", 2.315e-11)


register(
    [
        jaskowiec_sukumar_02,
        jaskowiec_sukumar_03,
        jaskowiec_sukumar_04,
        jaskowiec_sukumar_05,
        jaskowiec_sukumar_06,
        jaskowiec_sukumar_07,
        jaskowiec_sukumar_08,
        jaskowiec_sukumar_09,
        jaskowiec_sukumar_10,
        jaskowiec_sukumar_11,
        jaskowiec_sukumar_12,
        jaskowiec_sukumar_13,
        jaskowiec_sukumar_14,
        jaskowiec_sukumar_15,
        jaskowiec_sukumar_16,
        jaskowiec_sukumar_17,
        jaskowiec_sukumar_18,
        jaskowiec_sukumar_19a,
        jaskowiec_sukumar_19b,
        jaskowiec_sukumar_20,
    ]
)
