import json
import os

from ...helpers import article
from .._helpers import T2Scheme, untangle2

citation = article(
    authors=["F.D. Witherden", "P.E. Vincent"],
    title="On the identification of symmetric quadrature rules for finite element methods",
    journal="Computers & Mathematics with Applications",
    volume="69",
    number="10",
    month="may",
    year="2015",
    pages="1232–1241",
    url="https://doi.org/10.1016/j.camwa.2015.03.017",
)


def _read(filename):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    points, weights = untangle2(data)
    return weights, points, degree, citation


def witherden_vincent_01():
    return T2Scheme("Witherden-Vincent 1", *_read("wv01.json"))


def witherden_vincent_02():
    return T2Scheme("Witherden-Vincent 2", *_read("wv02.json"))


def witherden_vincent_04():
    return T2Scheme("Witherden-Vincent 4", *_read("wv04.json"))


def witherden_vincent_05():
    return T2Scheme("Witherden-Vincent 5", *_read("wv05.json"))


def witherden_vincent_06():
    return T2Scheme("Witherden-Vincent 6", *_read("wv06.json"))


def witherden_vincent_07():
    return T2Scheme("Witherden-Vincent 7", *_read("wv07.json"))


def witherden_vincent_08():
    return T2Scheme("Witherden-Vincent 8", *_read("wv08.json"))


def witherden_vincent_09():
    return T2Scheme("Witherden-Vincent 9", *_read("wv09.json"))


def witherden_vincent_10():
    return T2Scheme("Witherden-Vincent 10", *_read("wv10.json"))


def witherden_vincent_11():
    return T2Scheme("Witherden-Vincent 11", *_read("wv11.json"))


def witherden_vincent_12():
    return T2Scheme("Witherden-Vincent 12", *_read("wv12.json"))


def witherden_vincent_13():
    return T2Scheme("Witherden-Vincent 13", *_read("wv13.json"))


def witherden_vincent_14():
    return T2Scheme("Witherden-Vincent 14", *_read("wv14.json"))


def witherden_vincent_15():
    return T2Scheme("Witherden-Vincent 15", *_read("wv15.json"))


def witherden_vincent_16():
    return T2Scheme("Witherden-Vincent 16", *_read("wv16.json"))


def witherden_vincent_17():
    return T2Scheme("Witherden-Vincent 17", *_read("wv17.json"))


def witherden_vincent_18():
    return T2Scheme("Witherden-Vincent 18", *_read("wv18.json"))


def witherden_vincent_19():
    return T2Scheme("Witherden-Vincent 19", *_read("wv19.json"))


def witherden_vincent_20():
    return T2Scheme("Witherden-Vincent 20", *_read("wv20.json"))
