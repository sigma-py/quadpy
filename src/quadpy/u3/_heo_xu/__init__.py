import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["Sangwoo Heo", "Yuan Xu"],
    title="Constructing Fully Symmetric Cubature Formulae for the Sphere",
    journal="Mathematics of Computation",
    volume="70",
    number="233",
    month="jan",
    year="2001",
    pages="269-279",
    url="https://doi.org/10.1090/S0025-5718-00-01198-4",
)

this_dir = pathlib.Path(__file__).resolve().parent


def heo_xu_13():
    return _read(this_dir / "heo_xu_13.json", source)


def heo_xu_15():
    return _read(this_dir / "heo_xu_15.json", source)


def heo_xu_17():
    return _read(this_dir / "heo_xu_17.json", source)


def heo_xu_19a():
    return _read(this_dir / "heo_xu_19a.json", source)


def heo_xu_19b():
    return _read(this_dir / "heo_xu_19b.json", source)


def heo_xu_21a():
    return _read(this_dir / "heo_xu_21a.json", source)


def heo_xu_21b():
    return _read(this_dir / "heo_xu_21b.json", source)


def heo_xu_21c():
    return _read(this_dir / "heo_xu_21c.json", source)


def heo_xu_21d():
    return _read(this_dir / "heo_xu_21d.json", source)


def heo_xu_21e():
    return _read(this_dir / "heo_xu_21e.json", source)


def heo_xu_21f():
    return _read(this_dir / "heo_xu_21f.json", source)


def heo_xu_23a():
    return _read(this_dir / "heo_xu_23a.json", source)


def heo_xu_23b():
    return _read(this_dir / "heo_xu_23b.json", source)


def heo_xu_23c():
    return _read(this_dir / "heo_xu_23c.json", source)


def heo_xu_25a():
    return _read(this_dir / "heo_xu_25a.json", source)


def heo_xu_25b():
    return _read(this_dir / "heo_xu_25b.json", source)


def heo_xu_27a():
    return _read(this_dir / "heo_xu_27a.json", source)


def heo_xu_27b():
    return _read(this_dir / "heo_xu_27b.json", source)


def heo_xu_27c():
    return _read(this_dir / "heo_xu_27c.json", source)


def heo_xu_29():
    return _read(this_dir / "heo_xu_29.json", source)


def heo_xu_31():
    return _read(this_dir / "heo_xu_31.json", source)


def heo_xu_33():
    return _read(this_dir / "heo_xu_33.json", source)


def heo_xu_35():
    return _read(this_dir / "heo_xu_35.json", source)


def heo_xu_37():
    return _read(this_dir / "heo_xu_37.json", source)


def heo_xu_39a():
    return _read(this_dir / "heo_xu_39a.json", source)


def heo_xu_39b():
    return _read(this_dir / "heo_xu_39b.json", source)


register(
    [
        heo_xu_13,
        heo_xu_15,
        heo_xu_17,
        heo_xu_19a,
        heo_xu_19b,
        heo_xu_21a,
        heo_xu_21b,
        heo_xu_21c,
        heo_xu_21d,
        heo_xu_21e,
        heo_xu_21f,
        heo_xu_23a,
        heo_xu_23b,
        heo_xu_23c,
        heo_xu_25a,
        heo_xu_25b,
        heo_xu_27a,
        heo_xu_27b,
        heo_xu_27c,
        heo_xu_29,
        heo_xu_31,
        heo_xu_33,
        heo_xu_35,
        heo_xu_37,
        heo_xu_39a,
        heo_xu_39b,
    ]
)
