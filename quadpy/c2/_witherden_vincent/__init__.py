"""
"""
import json
import os

from ...helpers import article
from .._helpers import QuadrilateralScheme, concat, symm_r0, symm_s, symm_s_t, zero

_citation = article(
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


def _read(degree):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"wv{degree:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    assert degree == data.pop("degree")

    d = []
    if "zero" in data:
        d += [zero(data["zero"][0][0])]
    if "symm_r0" in data:
        d += [symm_r0(*data["symm_r0"])]
    if "symm_s" in data:
        d += [symm_s(*data["symm_s"])]
    if "symm_s_t" in data:
        d += [symm_s_t(*data["symm_s_t"])]

    weights, points = concat(*d)
    return QuadrilateralScheme(
        f"Witherden-Vincent {degree}", weights, points, degree, _citation
    )


def witherden_vincent_01():
    return _read(1)


def witherden_vincent_03():
    return _read(3)


def witherden_vincent_05():
    return _read(5)


def witherden_vincent_07():
    # TODO find error
    return _read(7)


def witherden_vincent_09():
    return _read(9)


def witherden_vincent_11():
    return _read(11)


def witherden_vincent_13():
    return _read(13)


def witherden_vincent_15():
    return _read(15)


def witherden_vincent_17():
    return _read(17)


def witherden_vincent_19():
    return _read(19)


def witherden_vincent_21():
    return _read(21)
