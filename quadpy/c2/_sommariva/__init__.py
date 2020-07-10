import json
import os

import numpy

from ...helpers import online
from .._helpers import C2Scheme

_source = online(
    authors=["Alvise Sommariva"],
    year="2012",
    url="http://www.math.unipd.it/~alvise/sets.html",
)


def _read(index, tol=1.0e-14):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"sommariva_{index:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    data = numpy.array(data["data"])
    points = data[:, :2]
    weights = data[:, 2]
    weights /= 4
    return C2Scheme(f"Sommariva {index}", weights, points, degree, _source, tol)


def sommariva_01():
    return _read(1)


def sommariva_02():
    return _read(2)


def sommariva_03():
    return _read(3)


def sommariva_04():
    return _read(4)


def sommariva_05():
    return _read(5)


def sommariva_06():
    return _read(6)


def sommariva_07():
    return _read(7)


def sommariva_08():
    return _read(8)


def sommariva_09():
    return _read(9)


def sommariva_10():
    return _read(10)


def sommariva_11():
    return _read(11)


def sommariva_12():
    return _read(12)


def sommariva_13():
    return _read(13)


def sommariva_14():
    return _read(14)


def sommariva_15():
    return _read(15)


def sommariva_16():
    return _read(16)


def sommariva_17():
    return _read(17)


def sommariva_18():
    return _read(18)


def sommariva_19():
    return _read(19)


def sommariva_20():
    return _read(20)


def sommariva_21():
    return _read(21, 1.922e-14)


def sommariva_22():
    return _read(22)


def sommariva_23():
    return _read(23)


def sommariva_24():
    return _read(24)


def sommariva_25():
    return _read(25)


def sommariva_26():
    return _read(26)


def sommariva_27():
    return _read(27)


def sommariva_28():
    return _read(28)


def sommariva_29():
    return _read(29)


def sommariva_30():
    return _read(30)


def sommariva_31():
    return _read(31)


def sommariva_32():
    return _read(32)


def sommariva_33():
    return _read(33)


def sommariva_34():
    return _read(34, 1.222e-14)


def sommariva_35():
    return _read(35)


def sommariva_36():
    return _read(36)


def sommariva_37():
    return _read(37, 1.306e-14)


def sommariva_38():
    return _read(38)


def sommariva_39():
    return _read(39)


def sommariva_40():
    return _read(40)


def sommariva_41():
    return _read(41)


def sommariva_42():
    return _read(42)


def sommariva_43():
    return _read(43)


def sommariva_44():
    return _read(44)


def sommariva_45():
    return _read(45)


def sommariva_46():
    return _read(46)


def sommariva_47():
    return _read(47)


def sommariva_48():
    return _read(48)


def sommariva_49():
    return _read(49)


def sommariva_50():
    return _read(50)


def sommariva_51():
    return _read(51)


def sommariva_52():
    return _read(52)


def sommariva_53():
    return _read(53)


def sommariva_54():
    return _read(54)


def sommariva_55():
    return _read(55)
