# -*- coding: utf-8 -*-
#
from .beckers_haegemans import BeckersHaegemans
from .gatermann import Gatermann
from .grundmann_moeller import GrundmannMoeller
from .hammer_marlowe_stroud import HammerMarloweStroud
from .hammer_stroud import HammerStroud
from .keast import Keast
from .liu_vinokur import LiuVinokur
from .maeztu_sainz import MaeztuSainz
from .newton_cotes import NewtonCotesClosed, NewtonCotesOpen
from .stroud import Stroud
from .shunn_ham import ShunnHam
from .vioreanu_rokhlin import VioreanuRokhlin
from .xiao_gimbutas import XiaoGimbutas
from .yu import Yu
from .zhang_cui_liu import ZhangCuiLiu
from .zienkiewicz import Zienkiewicz
from .walkington import Walkington
from .williams_shunn_jameson import WilliamsShunnJameson
from .witherden_vincent import WitherdenVincent

from .tools import show

from ..nsimplex import integrate

__all__ = [
    "BeckersHaegemans",
    "Gatermann",
    "GrundmannMoeller",
    "HammerMarloweStroud",
    "HammerStroud",
    "Keast",
    "LiuVinokur",
    "MaeztuSainz",
    "NewtonCotesClosed",
    "NewtonCotesOpen",
    "Stroud",
    "ShunnHam",
    "VioreanuRokhlin",
    "XiaoGimbutas",
    "Yu",
    "ZhangCuiLiu",
    "Zienkiewicz",
    "Walkington",
    "WilliamsShunnJameson",
    "WitherdenVincent",
    "show",
    "integrate",
]
