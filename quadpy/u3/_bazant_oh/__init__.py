import json
import pathlib

from ...helpers import article
from .._helpers import U3Scheme, cartesian_to_spherical, untangle2

source = article(
    authors=["P. Bažant", "B.H. Oh"],
    title="Efficient Numerical Integration on the Surface of a Sphere",
    journal="ZAMM",
    volume="66",
    number="1",
    year="1986",
    pages="37-49",
    url="https://doi.org/10.1002/zamm.19860660108",
)


def _read(filename):
    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / filename, "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    name = data.pop("name")
    tol = data.pop("test_tolerance")

    points, weights = untangle2(data)
    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, tol=tol)


def bazant_oh_09():
    return _read("bazant_oh_009.json")


def bazant_oh_11():
    return _read("bazant_oh_011.json")


def bazant_oh_13():
    return _read("bazant_oh_013.json")
