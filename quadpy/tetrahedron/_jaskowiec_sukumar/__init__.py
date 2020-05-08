import json
import os

import numpy

from ...helpers import article
from .._helpers import TetrahedronScheme

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


def _read(degree):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"js{degree:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    points = numpy.array(data["points"])
    points = numpy.column_stack([points, 1.0 - numpy.sum(points, axis=1)])
    weights = numpy.array(data["weights"])

    return TetrahedronScheme(
        f"Jaśkowiec-Sukumar {degree}", weights, points, degree, citation
    )


def jaskowiec_sukumar_02():
    return _read(2)
