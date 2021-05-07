# ENH quadpy-optimize
import pathlib

from ...helpers import article
from .._helpers import _read, register

_source = article(
    authors=["KyoungJoong Kim", "ManSuk Song"],
    title="Symmetric quadrature formulas over a unit disk",
    journal="Korean J. Comp. & Appl. Math.",
    year="1997",
    volume="4",
    pages="179-192",
    url="https://doi.org/10.1007/BF03011388",
)

this_dir = pathlib.Path(__file__).resolve().parent


def kim_song_1():
    return _read(this_dir / "kim_song_01.json", _source)


def kim_song_2():
    return _read(this_dir / "kim_song_02.json", _source)


def kim_song_3():
    return _read(this_dir / "kim_song_03.json", _source)


def kim_song_4():
    return _read(this_dir / "kim_song_04.json", _source)


def kim_song_5():
    return _read(this_dir / "kim_song_05.json", _source)


def kim_song_6():
    return _read(this_dir / "kim_song_06.json", _source)


def kim_song_7():
    return _read(this_dir / "kim_song_07.json", _source)


# TODO find issue
def kim_song_8():
    return _read(this_dir / "kim_song_08.json", _source)


def kim_song_9():
    return _read(this_dir / "kim_song_09.json", _source)


def kim_song_10():
    return _read(this_dir / "kim_song_10.json", _source)


def kim_song_11():
    return _read(this_dir / "kim_song_11.json", _source)


def kim_song_12():
    return _read(this_dir / "kim_song_12.json", _source)


def kim_song_13():
    return _read(this_dir / "kim_song_13.json", _source)


def kim_song_14():
    return _read(this_dir / "kim_song_14.json", _source)


def kim_song_15():
    return _read(this_dir / "kim_song_15.json", _source)


register(
    [
        kim_song_1,
        kim_song_2,
        kim_song_3,
        kim_song_4,
        kim_song_5,
        kim_song_6,
        kim_song_7,
        kim_song_8,
        kim_song_9,
        kim_song_10,
        kim_song_11,
        kim_song_12,
        kim_song_13,
        kim_song_14,
        kim_song_15,
    ]
)
