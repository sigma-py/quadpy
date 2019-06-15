# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from .stroud_secrest import StroudSecrest


Stroud = {
    "5-1": StroudSecrest["VII"],
    "5-2": StroudSecrest["VIII"],
    "5-3": StroudSecrest["IX"],
    # '7-1': StroudSecrest["X"],
    "7-2": StroudSecrest["XI"],
}
