import json
import os
import re

from ...helpers import article
from .._helpers import SphereScheme, cartesian_to_spherical, untangle2

citation = article(
    authors=["P. Bažant", "B.H. Oh"],
    title="Efficient Numerical Integration on the Surface of a Sphere",
    journal="ZAMM",
    volume="66",
    number="1",
    year="1986",
    pages="37-49",
    url="https://doi.org/10.1002/zamm.19860660108",
)


def _read(index):
    name = f"BazantOh({index})"

    this_dir = os.path.dirname(os.path.realpath(__file__))

    m = re.match("([0-9]+)([a-z]*)", index)
    filename = "bazant_oh_{:03d}{}.json".format(int(m.group(1)), m.group(2))
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    points, weights = untangle2(data)
    azimuthal_polar = cartesian_to_spherical(points)
    return SphereScheme(name, weights, points, azimuthal_polar, degree, citation)


def bazant_oh_09():
    return _read("9")


def bazant_oh_11():
    return _read("11")


def bazant_oh_13():
    return _read("13")
