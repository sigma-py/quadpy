import json
import os

import numpy

from ...helpers import article
from .._helpers import T3Scheme

citation = article(
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


def _read(string):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"js{string}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    points = numpy.array(data["points"])
    points = numpy.column_stack([points, 1.0 - numpy.sum(points, axis=1)])
    weights = numpy.array(data["weights"])

    return T3Scheme(f"Jaśkowiec-Sukumar {string}", weights, points, degree, citation)


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
    return _read("15")


def jaskowiec_sukumar_16():
    return _read("16")


def jaskowiec_sukumar_17():
    return _read("17")


def jaskowiec_sukumar_18():
    return _read("18")


def jaskowiec_sukumar_19a():
    return _read("19a")


def jaskowiec_sukumar_19b():
    return _read("19b")


def jaskowiec_sukumar_20():
    return _read("20")
