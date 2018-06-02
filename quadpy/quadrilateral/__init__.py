# -*- coding: utf-8 -*-
#

from .albrecht_collatz import AlbrechtCollatz
from .burnside import Burnside
from .cools_haegemans_1985 import CoolsHaegemans1985
from .cools_haegemans_1988 import CoolsHaegemans1988
from .dunavant import Dunavant
from .hammer_stroud import HammerStroud
from .irwin import Irwin
from .maxwell import Maxwell
from .meister import Meister
from .miller import Miller
from .morrow_patterson import MorrowPatterson
from .phillips import Phillips
from .rabinowitz_richter import RabinowitzRichter
from .stroud import Stroud
from .tyler import Tyler
from .wissmann_becker import WissmannBecker

from .product import Product
from .stroudn import StroudN

from .tools import show, plot

__all__ = [
    "AlbrechtCollatz",
    "Burnside",
    "CoolsHaegemans1985",
    "CoolsHaegemans1988",
    "Dunavant",
    "HammerStroud",
    "Irwin",
    "Maxwell",
    "Meister",
    "Miller",
    "MorrowPatterson",
    "Phillips",
    "RabinowitzRichter",
    "Stroud",
    "Tyler",
    "WissmannBecker",
    "Product",
    "StroudN",
    "show",
    "plot",
]
