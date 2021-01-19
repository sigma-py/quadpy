import pathlib

import numpy as np
from sympy import Rational as frac

from ...helpers import article, book
from .._hammer_marlowe_stroud import hammer_marlowe_stroud_3 as strang_fix_cowper_02
from .._helpers import T2Scheme, _read, register

source = book(
    authors=["Gilbert Strang", "George Fix"],
    title="An Analysis of the Finite Element Method",
    publisher="Wellesley-Cambridge Press",
    year="1973",
    isbn="096140888X",
    url="https://bookstore.siam.org/wc08/",
)

source = article(
    authors=["G.R. Cowper"],
    title="Gaussian quadrature formulas for triangles",
    journal="Numerical Methods in Engineering",
    volume="7",
    number="3",
    year="1973",
    pages="405–408",
    url="https://doi.org/10.1002/nme.1620070316",
)
# See also
# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html

this_dir = pathlib.Path(__file__).resolve().parent


def strang_fix_cowper_01():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Strang-Fix-Cowper 1", d, 2, source)


def strang_fix_cowper_03():
    d = {"centroid": [[-frac(9, 16)]], "d3_aa": [[frac(25, 48)], [frac(1, 5)]]}
    return T2Scheme("Strang-Fix-Cowper 3", d, 3, source)


def strang_fix_cowper_04():
    roots = np.polynomial.polynomial.polyroots([-1, 15, -60, 60])
    d = {"d3_ab": [[1 / 6], [roots[2]], [roots[1]]]}
    return T2Scheme("Strang-Fix-Cowper 4", d, 3, source, 2.463e-15)


def strang_fix_cowper_05():
    return _read(this_dir / "strang_fix_cowper_05.json", source)


def strang_fix_cowper_06():
    d = {
        "centroid": [[3 / 8]],
        "d3_ab": [[5 / 48], [0.736712498968435], [0.237932366472434]],
    }
    return T2Scheme("Strang-Fix-Cowper 6", d, 4, source, 5.440e-15)


def strang_fix_cowper_07():
    d = {
        "centroid": [[9 / 40]],
        "d3_aa": [
            [0.12593918054482717, 0.13239415278850616],
            [0.10128650732345633, 0.47014206410511505],
        ],
    }
    return T2Scheme("Strang-Fix-Cowper 7", d, 5, source, 5.551e-16)


def strang_fix_cowper_08():
    return _read(this_dir / "strang_fix_cowper_08.json", source)


def strang_fix_cowper_09():
    return _read(this_dir / "strang_fix_cowper_09.json", source)


def strang_fix_cowper_10():
    return _read(this_dir / "strang_fix_cowper_10.json", source)


register(
    [
        strang_fix_cowper_01,
        strang_fix_cowper_02,
        strang_fix_cowper_03,
        strang_fix_cowper_04,
        strang_fix_cowper_05,
        strang_fix_cowper_06,
        strang_fix_cowper_07,
        strang_fix_cowper_08,
        strang_fix_cowper_09,
        strang_fix_cowper_10,
    ]
)
