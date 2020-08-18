import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["B. Vioreanu", "V. Rokhlin"],
    title="Spectra of Multiplication Operators as a Numerical Tool",
    year="2014",
    journal="SIAM J. Sci. Comput.",
    volume="36",
    number="1",
    url="https://doi.org/10.1137/110860082",
)
# http://www.cs.yale.edu/publications/techreports/tr1443.pdf

# Data adapted from modepy
# https://github.com/inducer/modepy/blob/master/modepy/quadrature/vr_quad_data_tri.py

this_dir = pathlib.Path(__file__).resolve().parent


def vioreanu_rokhlin_0():
    return _read(this_dir / "vr00.json", source)


def vioreanu_rokhlin_1():
    return _read(this_dir / "vr01.json", source)


def vioreanu_rokhlin_2():
    return _read(this_dir / "vr02.json", source)


def vioreanu_rokhlin_3():
    return _read(this_dir / "vr03.json", source)


def vioreanu_rokhlin_4():
    return _read(this_dir / "vr04.json", source)


def vioreanu_rokhlin_5():
    return _read(this_dir / "vr05.json", source)


def vioreanu_rokhlin_6():
    return _read(this_dir / "vr06.json", source)


def vioreanu_rokhlin_7():
    return _read(this_dir / "vr07.json", source)


def vioreanu_rokhlin_8():
    return _read(this_dir / "vr08.json", source)


def vioreanu_rokhlin_9():
    return _read(this_dir / "vr09.json", source)


register(
    [
        vioreanu_rokhlin_0,
        vioreanu_rokhlin_1,
        vioreanu_rokhlin_2,
        vioreanu_rokhlin_3,
        vioreanu_rokhlin_4,
        vioreanu_rokhlin_5,
        vioreanu_rokhlin_6,
        vioreanu_rokhlin_7,
        vioreanu_rokhlin_8,
        vioreanu_rokhlin_9,
    ]
)
