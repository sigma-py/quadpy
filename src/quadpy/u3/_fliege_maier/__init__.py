import pathlib

from ...helpers import online
from .._helpers import _read, register

source = online(
    authors=["Jörg Fliege", "Ulrike Maier"],
    title="A Two-Stage Approach for Computing Cubature Formulae for the Sphere",
    url="https://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html",
)

this_dir = pathlib.Path(__file__).resolve().parent


def fliege_maier_04():
    return _read(this_dir / "fliege_maier_04.json", source)


def fliege_maier_09():
    return _read(this_dir / "fliege_maier_09.json", source)


def fliege_maier_16():
    return _read(this_dir / "fliege_maier_16.json", source)


def fliege_maier_25():
    return _read(this_dir / "fliege_maier_25.json", source)


register([fliege_maier_04, fliege_maier_09, fliege_maier_16, fliege_maier_25])
