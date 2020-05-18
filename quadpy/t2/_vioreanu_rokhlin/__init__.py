import json
import os

from ...helpers import article
from .._helpers import T2Scheme, untangle2

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


def _read(index, tol=1.0e-14):
    # Use JSON instead of YAML because
    #  * a parse is in the Python standard library, and
    #  * it's _much_ faster to parse <https://stackoverflow.com/a/50685946/353337>
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"vr{index:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    points, weights = untangle2(data)
    weights /= 2
    return T2Scheme(f"Vioreanu-Rokhlin {index}", weights, points, degree, source, tol)


def vioreanu_rokhlin_00():
    return _read(0)


def vioreanu_rokhlin_01():
    return _read(1)


def vioreanu_rokhlin_02():
    return _read(2)


def vioreanu_rokhlin_03():
    return _read(3)


def vioreanu_rokhlin_04():
    return _read(4)


def vioreanu_rokhlin_05():
    return _read(5)


def vioreanu_rokhlin_06():
    return _read(6)


def vioreanu_rokhlin_07():
    return _read(7)


def vioreanu_rokhlin_08():
    return _read(8)


def vioreanu_rokhlin_09():
    return _read(9)


def vioreanu_rokhlin_10():
    return _read(10)


def vioreanu_rokhlin_11():
    return _read(11)


def vioreanu_rokhlin_12():
    return _read(12)


def vioreanu_rokhlin_13():
    return _read(13, 3.045e-14)


def vioreanu_rokhlin_14():
    return _read(14, 2.055e-14)


def vioreanu_rokhlin_15():
    return _read(15, 2.112e-14)


def vioreanu_rokhlin_16():
    return _read(16, 2.523e-13)


def vioreanu_rokhlin_17():
    return _read(17)


def vioreanu_rokhlin_18():
    return _read(18, 2.895e-14)


def vioreanu_rokhlin_19():
    return _read(19, 4.085e-12)
