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


def vioreanu_rokhlin_00():
    return _read(this_dir / "vr00.json", source)


def vioreanu_rokhlin_01():
    return _read(this_dir / "vr01.json", source)


def vioreanu_rokhlin_02():
    return _read(this_dir / "vr02.json", source)


def vioreanu_rokhlin_03():
    return _read(this_dir / "vr03.json", source)


def vioreanu_rokhlin_04():
    return _read(this_dir / "vr04.json", source)


def vioreanu_rokhlin_05():
    return _read(this_dir / "vr05.json", source)


def vioreanu_rokhlin_06():
    return _read(this_dir / "vr06.json", source)


def vioreanu_rokhlin_07():
    return _read(this_dir / "vr07.json", source)


def vioreanu_rokhlin_08():
    return _read(this_dir / "vr08.json", source)


def vioreanu_rokhlin_09():
    return _read(this_dir / "vr09.json", source)


def vioreanu_rokhlin_10():
    return _read(this_dir / "vr10.json", source)


def vioreanu_rokhlin_11():
    return _read(this_dir / "vr11.json", source)


def vioreanu_rokhlin_12():
    return _read(this_dir / "vr12.json", source)


def vioreanu_rokhlin_13():
    return _read(this_dir / "vr13.json", source)


def vioreanu_rokhlin_14():
    return _read(this_dir / "vr14.json", source)


def vioreanu_rokhlin_15():
    return _read(this_dir / "vr15.json", source)


def vioreanu_rokhlin_16():
    return _read(this_dir / "vr16.json", source)


def vioreanu_rokhlin_17():
    return _read(this_dir / "vr17.json", source)


def vioreanu_rokhlin_18():
    return _read(this_dir / "vr18.json", source)


def vioreanu_rokhlin_19():
    return _read(this_dir / "vr19.json", source)


register(
    [
        vioreanu_rokhlin_00,
        vioreanu_rokhlin_01,
        vioreanu_rokhlin_02,
        vioreanu_rokhlin_03,
        vioreanu_rokhlin_04,
        vioreanu_rokhlin_05,
        vioreanu_rokhlin_06,
        vioreanu_rokhlin_07,
        vioreanu_rokhlin_08,
        vioreanu_rokhlin_09,
        vioreanu_rokhlin_10,
        vioreanu_rokhlin_10,
        vioreanu_rokhlin_11,
        vioreanu_rokhlin_12,
        vioreanu_rokhlin_13,
        vioreanu_rokhlin_14,
        vioreanu_rokhlin_15,
        vioreanu_rokhlin_16,
        vioreanu_rokhlin_17,
        vioreanu_rokhlin_18,
        vioreanu_rokhlin_19,
    ]
)
