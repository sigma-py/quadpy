# -*- coding: utf-8 -*-
#
from ._grundmann_moeller import grundmann_moeller
from ._hammer_stroud import hammer_stroud_1a, hammer_stroud_1b, hammer_stroud_2
from ._stroud import (
    stroud_tn_1_1,
    stroud_tn_1_2,
    stroud_tn_2_1a,
    stroud_tn_2_1b,
    stroud_tn_2_2,
    stroud_tn_3_1,
    stroud_tn_3_2,
    stroud_tn_3_3,
    stroud_tn_3_4,
    stroud_tn_3_5,
    stroud_tn_3_6a,
    stroud_tn_3_6b,
    stroud_tn_3_7,
    stroud_tn_3_8,
    stroud_tn_3_9,
    stroud_tn_3_10,
    stroud_tn_3_11,
    stroud_tn_4_1,
    stroud_tn_5_1,
    stroud_tn_5_2,
)
from ._stroud_1961 import stroud_1961
from ._stroud_1964 import stroud_1964a, stroud_1964b
from ._stroud_1966 import (
    stroud_1966_i,
    stroud_1966_ii,
    stroud_1966_iii,
    stroud_1966_iv,
    stroud_1966_v,
    stroud_1966_vi,
    stroud_1966_vii,
)
from ._stroud_1969 import stroud_1969
from ._walkington import (
    walkington_1,
    walkington_2,
    walkington_3,
    walkington_5,
    walkington_7,
)

from ._helpers import transform, get_vol, NSimplexScheme

__all__ = [
    "grundmann_moeller",
    "hammer_stroud_1a",
    "hammer_stroud_1b",
    "hammer_stroud_2",
    "stroud_1961",
    "stroud_1964a",
    "stroud_1964b",
    "stroud_1966_i",
    "stroud_1966_ii",
    "stroud_1966_iii",
    "stroud_1966_iv",
    "stroud_1966_v",
    "stroud_1966_vi",
    "stroud_1966_vii",
    "stroud_1969",
    "stroud_tn_1_1",
    "stroud_tn_1_2",
    "stroud_tn_2_1a",
    "stroud_tn_2_1b",
    "stroud_tn_2_2",
    "stroud_tn_3_1",
    "stroud_tn_3_2",
    "stroud_tn_3_3",
    "stroud_tn_3_4",
    "stroud_tn_3_5",
    "stroud_tn_3_6a",
    "stroud_tn_3_6b",
    "stroud_tn_3_7",
    "stroud_tn_3_8",
    "stroud_tn_3_9",
    "stroud_tn_3_10",
    "stroud_tn_3_11",
    "stroud_tn_4_1",
    "stroud_tn_5_1",
    "stroud_tn_5_2",
    "walkington_1",
    "walkington_2",
    "walkington_3",
    "walkington_5",
    "walkington_7",
    #
    "transform",
    "get_vol",
    "NSimplexScheme",
]
