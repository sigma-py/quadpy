import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["Johannes W. Wissmann", "Thomas Becker"],
    title="Partially Symmetric Cubature Formulas for Even Degrees of Exactness",
    journal="SIAM J. Numer. Anal.",
    volume="23",
    number="3",
    pages="676–685",
    url="https://doi.org/10.1137/0723043",
)

this_dir = pathlib.Path(__file__).resolve().parent


def wissmann_becker_4_1():
    return _read(this_dir / "wissmann_becker_4_1.json", source)


def wissmann_becker_4_2():
    return _read(this_dir / "wissmann_becker_4_2.json", source)


def wissmann_becker_6_1():
    return _read(this_dir / "wissmann_becker_6_1.json", source)


def wissmann_becker_6_2():
    return _read(this_dir / "wissmann_becker_6_2.json", source)


def wissmann_becker_8_1():
    return _read(this_dir / "wissmann_becker_8_1.json", source)


def wissmann_becker_8_2():
    return _read(this_dir / "wissmann_becker_8_2.json", source)


register(
    [
        wissmann_becker_4_1,
        wissmann_becker_4_2,
        wissmann_becker_6_1,
        wissmann_becker_6_2,
        wissmann_becker_8_1,
        wissmann_becker_8_2,
    ]
)
