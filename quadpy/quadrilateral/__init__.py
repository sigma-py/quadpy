# -*- coding: utf-8 -*-
#

from .albrecht_collatz import AlbrechtCollatz
from .burnside import Burnside
from .cohen_gismalla import CohenGismalla
from .cools_haegemans_1985 import CoolsHaegemans1985
from .cools_haegemans_1988 import CoolsHaegemans1988
from .dunavant import Dunavant
from .hammer_stroud import HammerStroud
from .haegemans_piessens import HaegemansPiessens
from .irwin import Irwin
from .maxwell import Maxwell
from .meister import Meister
from .miller import Miller
from .morrow_patterson import MorrowPatterson
from .phillips import Phillips
from .piessens_haegemans import PiessensHaegemans
from .rabinowitz_richter import RabinowitzRichter
from .schmid import Schmid
from .sommariva import Sommariva
from .stroud import Stroud
from .tyler import Tyler
from .waldron import Waldron
from .wissmann_becker import WissmannBecker
from .witherden_vincent import WitherdenVincent

from .product import Product
from .stroudn import StroudN

from .tools import show, plot

from ..ncube import transform, integrate
from ..ncube import ncube_points as rectangle_points

__all__ = [
    "AlbrechtCollatz",
    "Burnside",
    "CohenGismalla",
    "CoolsHaegemans1985",
    "CoolsHaegemans1988",
    "Dunavant",
    "HammerStroud",
    "HaegemansPiessens",
    "Irwin",
    "Maxwell",
    "Meister",
    "Miller",
    "MorrowPatterson",
    "PiessensHaegemans",
    "Phillips",
    "RabinowitzRichter",
    "Stroud",
    "Tyler",
    "Waldron",
    "WissmannBecker",
    "WitherdenVincent",
    "Product",
    "Schmid",
    "Sommariva",
    "StroudN",
    "show",
    "plot",
    "transform",
    "integrate",
    "rectangle_points",
]
