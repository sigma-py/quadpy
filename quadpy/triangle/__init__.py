# -*- coding: utf-8 -*-
#
from .albrecht_collatz import AlbrechtCollatz
from .centroid import Centroid
from .cools_haegemans import CoolsHaegemans
from .cubtri import Cubtri
from .berntsen_espelid import BerntsenEspelid
from .dunavant import Dunavant
from .gatermann import Gatermann
from .grundmann_moeller import GrundmannMoeller
from .hammer_marlowe_stroud import HammerMarloweStroud
from .hammer_stroud import HammerStroud
from .hillion import Hillion
from .laursen_gellert import LaursenGellert
from .lether import Lether
from .liu_vinokur import LiuVinokur
from .lyness_jespersen import LynessJespersen
from .newton_cotes import NewtonCotesClosed, NewtonCotesOpen
from .papanicolopulos import Papanicolopulos
from .rathod_nagaraja_venkatesudu import RathodNagarajaVenkatesudu
from .seven_point import SevenPoint
from .strang import Strang
from .stroud import Stroud
from .taylor_wingate_bos import TaylorWingateBos
from .triex import Triex
from .vertex import Vertex
from .vioreanu_rokhlin import VioreanuRokhlin
from .walkington import Walkington
from .wandzura_xiao import WandzuraXiao
from .williams_shunn_jameson import WilliamsShunnJameson
from .witherden_vincent import WitherdenVincent
from .xiao_gimbutas import XiaoGimbutas
from .zhang_cui_liu import ZhangCuiLiu

from .tools import show, plot, integrate_adaptive

from ..nsimplex import transform, get_vol, integrate

__all__ = [
    "AlbrechtCollatz",
    "Centroid",
    "CoolsHaegemans",
    "Cubtri",
    "BerntsenEspelid",
    "Dunavant",
    "Gatermann",
    "GrundmannMoeller",
    "HammerMarloweStroud",
    "HammerStroud",
    "Hillion",
    "LaursenGellert",
    "Lether",
    "LiuVinokur",
    "LynessJespersen",
    "NewtonCotesClosed",
    "NewtonCotesOpen",
    "Papanicolopulos",
    "RathodNagarajaVenkatesudu",
    "SevenPoint",
    "Strang",
    "Stroud",
    "TaylorWingateBos",
    "Triex",
    "Vertex",
    "VioreanuRokhlin",
    "Walkington",
    "WandzuraXiao",
    "WilliamsShunnJameson",
    "WitherdenVincent",
    "XiaoGimbutas",
    "ZhangCuiLiu",
    "show",
    "plot",
    "integrate_adaptive",
    "transform",
    "get_vol",
    "integrate",
]
