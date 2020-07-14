import pathlib

import numpy

from ...helpers import online
from .._helpers import _read

source = online(
    authors=["Jörg Fliege", "Ulrike Maier"],
    title="A Two-Stage Approach for Computing Cubature Formulae for the Sphere",
    url="https://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html",
)

this_dir = pathlib.Path(__file__).resolve().parent
weight_factor = 1 / 4 / numpy.pi


def fliege_maier_04():
    return _read(this_dir / "fliege_maier_004.json", source, weight_factor)


def fliege_maier_09():
    return _read(this_dir / "fliege_maier_009.json", source, weight_factor)


def fliege_maier_16():
    return _read(this_dir / "fliege_maier_016.json", source, weight_factor)


def fliege_maier_25():
    return _read(this_dir / "fliege_maier_025.json", source, weight_factor)
