# -*- coding: utf-8 -*-
#
from .albrecht_collatz import AlbrechtCollatz
from .bazant_oh import BazantOh
from .fliege_maier import FliegeMaier
from .heo_xu import HeoXu
from .lebedev import Lebedev
from .mclaren import McLaren
from .stroud import Stroud

from .tools import area, show, plot, integrate, integrate_spherical

__all__ = [
    "AlbrechtCollatz",
    "BazantOh",
    "FliegeMaier",
    "HeoXu",
    "Lebedev",
    "McLaren",
    "Stroud",
    "area",
    "show",
    "plot",
    "integrate",
    "integrate_spherical",
]
