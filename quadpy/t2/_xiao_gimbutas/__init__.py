import json
import os

from ...helpers import article
from .._helpers import T2Scheme, untangle2

source = article(
    authors=["Hong Xiao", "Zydrunas Gimbutas"],
    title="A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions",
    journal="Computers & Mathematics with Applications",
    volume="59",
    number="2",
    month="jan",
    year="2010",
    pages="663–676",
    url="https://doi.org/10.1016/j.camwa.2009.10.027",
)
# Data extracted from
# https://people.sc.fsu.edu/~jburkardt/f_src/triangle_symq_rule/triangle_symq_rule.f90


# TODO mpmath
def _read(degree):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"xg{degree:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    name = f"Xiao-Gimbutas {degree}"
    points, weights = untangle2(data)
    return T2Scheme(name, weights, points, degree, source)


def xiao_gimbutas_01():
    return _read(1)


def xiao_gimbutas_02():
    return _read(2)


def xiao_gimbutas_03():
    return _read(3)


def xiao_gimbutas_04():
    return _read(4)


def xiao_gimbutas_05():
    return _read(5)


def xiao_gimbutas_06():
    return _read(6)


def xiao_gimbutas_07():
    return _read(7)


def xiao_gimbutas_08():
    return _read(8)


def xiao_gimbutas_09():
    return _read(9)


def xiao_gimbutas_10():
    return _read(10)


def xiao_gimbutas_11():
    return _read(11)


def xiao_gimbutas_12():
    return _read(12)


def xiao_gimbutas_13():
    return _read(13)


def xiao_gimbutas_14():
    return _read(14)


def xiao_gimbutas_15():
    return _read(15)


def xiao_gimbutas_16():
    return _read(16)


def xiao_gimbutas_17():
    return _read(17)


def xiao_gimbutas_18():
    return _read(18)


def xiao_gimbutas_19():
    return _read(19)


def xiao_gimbutas_20():
    return _read(20)


def xiao_gimbutas_21():
    return _read(21)


def xiao_gimbutas_22():
    return _read(22)


def xiao_gimbutas_23():
    return _read(23)


def xiao_gimbutas_24():
    return _read(24)


def xiao_gimbutas_25():
    return _read(25)


def xiao_gimbutas_26():
    return _read(26)


def xiao_gimbutas_27():
    return _read(27)


def xiao_gimbutas_28():
    return _read(28)


def xiao_gimbutas_29():
    return _read(29)


def xiao_gimbutas_30():
    return _read(30)


def xiao_gimbutas_31():
    return _read(31)


def xiao_gimbutas_32():
    return _read(32)


def xiao_gimbutas_33():
    return _read(33)


def xiao_gimbutas_34():
    return _read(34)


def xiao_gimbutas_35():
    return _read(35)


def xiao_gimbutas_36():
    return _read(36)


def xiao_gimbutas_37():
    return _read(37)


def xiao_gimbutas_38():
    return _read(38)


def xiao_gimbutas_39():
    return _read(39)


def xiao_gimbutas_40():
    return _read(40)


def xiao_gimbutas_41():
    return _read(41)


def xiao_gimbutas_42():
    return _read(42)


def xiao_gimbutas_43():
    return _read(43)


def xiao_gimbutas_44():
    return _read(44)


def xiao_gimbutas_45():
    return _read(45)


def xiao_gimbutas_46():
    return _read(46)


def xiao_gimbutas_47():
    return _read(47)


def xiao_gimbutas_48():
    return _read(48)


def xiao_gimbutas_49():
    return _read(49)


def xiao_gimbutas_50():
    return _read(50)
