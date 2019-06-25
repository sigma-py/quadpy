# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from ._stroud_secrest import (
    stroud_secrest_vii as stroud_e3r_5_1,
    stroud_secrest_viii as stroud_e3r_5_2,
    stroud_secrest_ix as stroud_e3r_5_3,
    # stroud_secrest_x as stroud_e3r_7_1,
    stroud_secrest_xi as stroud_e3r_7_2,
)

__all__ = [
    "stroud_e3r_5_1",
    "stroud_e3r_5_2",
    "stroud_e3r_5_3",
    # "stroud_e3r_7_1",
    "stroud_e3r_7_2",
]
