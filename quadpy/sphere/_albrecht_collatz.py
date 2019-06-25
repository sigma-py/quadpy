# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._helpers import cartesian_to_spherical_sympy, SphereScheme
from ..helpers import untangle, pm_array0, fsd, pm, article

citation = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)


def albrecht_collatz_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r, s = [sqrt((5 + t * sqrt(5)) / 10) for t in [+1, -1]]
    data = [
        (frac(1, 12), pm_array0(3, [r, s], [0, 1])),
        (frac(1, 12), pm_array0(3, [r, s], [1, 2])),
        (frac(1, 12), pm_array0(3, [r, s], [2, 0])),
    ]

    points, weights = untangle(data)
    azimuthal_polar = cartesian_to_spherical_sympy(points)
    return SphereScheme(
        "Albrecht-Collatz 1", weights, points, azimuthal_polar, 5, citation
    )


def albrecht_collatz_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = 1
    s = sqrt(frac(1, 3))
    data = [(frac(8, 120), fsd(3, (r, 1))), (frac(9, 120), pm(3, s))]

    points, weights = untangle(data)
    azimuthal_polar = cartesian_to_spherical_sympy(points)
    return SphereScheme(
        "Albrecht-Collatz 2", weights, points, azimuthal_polar, 5, citation
    )


def albrecht_collatz_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = 1
    s = sqrt(frac(1, 2))
    data = [(frac(1, 30), fsd(3, (r, 1))), (frac(2, 30), fsd(3, (s, 2)))]

    points, weights = untangle(data)
    azimuthal_polar = cartesian_to_spherical_sympy(points)
    return SphereScheme(
        "Albrecht-Collatz 3", weights, points, azimuthal_polar, 5, citation
    )


def albrecht_collatz_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r, s = [sqrt((3 + t * sqrt(5)) / 6) for t in [+1, -1]]
    t = sqrt(frac(1, 3))
    data = [
        (frac(1, 20), pm_array0(3, [r, s], [0, 1])),
        (frac(1, 20), pm_array0(3, [r, s], [1, 2])),
        (frac(1, 20), pm_array0(3, [r, s], [2, 0])),
        (frac(1, 20), pm(3, t)),
    ]

    points, weights = untangle(data)
    azimuthal_polar = cartesian_to_spherical_sympy(points)
    return SphereScheme(
        "Albrecht-Collatz 4", weights, points, azimuthal_polar, 5, citation
    )


def albrecht_collatz_5(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = 1
    s = sqrt(frac(1, 2))
    t = sqrt(frac(1, 3))

    data = [
        (frac(40, 840), fsd(3, (r, 1))),
        (frac(32, 840), fsd(3, (s, 2))),
        (frac(27, 840), pm(3, t)),
    ]

    points, weights = untangle(data)
    azimuthal_polar = cartesian_to_spherical_sympy(points)
    return SphereScheme(
        "Albrecht-Collatz 5", weights, points, azimuthal_polar, 7, citation
    )
