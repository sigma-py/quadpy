import json
import os
import re
import warnings

import numpy

from ...helpers import online
from .._helpers import SphereScheme, cartesian_to_spherical

citation = online(
    authors=["Jörg Fliege", "Ulrike Maier"],
    title="A Two-Stage Approach for Computing Cubature Formulae for the Sphere",
    url="https://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html",
)


def _read(index):
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

    azimuthal_polar = cartesian_to_spherical(points)
    return SphereScheme(name, weights, points, azimuthal_polar, degree, citation)


def fliege_maier_04():
    return _read("4")


def fliege_maier_09():
    return _read("9")


def fliege_maier_16():
    return _read("16")


def fliege_maier_25():
    return _read("25")
