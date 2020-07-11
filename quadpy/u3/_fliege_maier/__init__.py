import json
import os
import re
import warnings

import numpy

from ...helpers import online
from .._helpers import U3Scheme, cartesian_to_spherical

source = online(
    authors=["Jörg Fliege", "Ulrike Maier"],
    title="A Two-Stage Approach for Computing Cubature Formulae for the Sphere",
    url="https://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html",
)


def _read(index, tol):
    warnings.warn("The Fliege-Maier schemes are only single-precision.")

    name = f"FliegeMaier({index})"

    this_dir = os.path.dirname(os.path.realpath(__file__))

    m = re.match("([0-9]+)([a-z]*)", index)
    filename = "fliege_maier_{:03d}{}.json".format(int(m.group(1)), m.group(2))
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    data = numpy.array(data["data"])
    points = data[:, :3]
    weights = data[:, 3] / 4 / numpy.pi

    theta_phi = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, theta_phi, degree, source, tol)


def fliege_maier_04():
    return _read("4", 1.564e-07)


def fliege_maier_09():
    return _read("9", 1.602e-12)


def fliege_maier_16():
    return _read("16", 7.773e-12)


def fliege_maier_25():
    return _read("25", 7.886e-12)
