# -*- coding: utf-8 -*-
#
from .albrecht_collatz import AlbrechtCollatz
from .heo_xu import HeoXu
from .lebedev import Lebedev
from .mclaren import McLaren
from .stroud import Stroud

from .tools import area, show, plot, integrate, integrate_spherical

__all__ = [
    "AlbrechtCollatz",
    "HeoXu ",
    "Lebedev",
    "McLaren",
    "Stroud",
    "area",
    "show",
    "plot",
    "integrate",
    "integrate_spherical",
]
