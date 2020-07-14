import json
import pathlib
import warnings

import numpy

from ...helpers import online
from .._helpers import U3Scheme, cartesian_to_spherical

source = online(
    authors=["Jörg Fliege", "Ulrike Maier"],
    title="A Two-Stage Approach for Computing Cubature Formulae for the Sphere",
    url="https://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html",
)


def _read(filename, tol):
    warnings.warn("The Fliege-Maier schemes are only single-precision.")

    this_dir = pathlib.Path(__file__).resolve().parent

    with open(this_dir / filename, "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    name = data.pop("name")
    tol = data.pop("test_tolerance")

    data = numpy.array(data["data"])
    points = data[:, :3]
    weights = data[:, 3] / 4 / numpy.pi

    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, tol)


def fliege_maier_04():
    return _read("fliege_maier_004.json")


def fliege_maier_09():
    return _read("fliege_maier_009.json")


def fliege_maier_16():
    return _read("fliege_maier_016.json")


def fliege_maier_25():
    return _read("fliege_maier_025.json")
