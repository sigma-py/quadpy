from . import _cools_haegemans
from . import _haegemans_piessens
from . import _rabinowitz_richter
from . import _stroud
from . import _stroud_secrest

from ._helpers import schemes, get_good_scheme

__all__ = [
    "_cools_haegemans",
    "_haegemans_piessens",
    "_rabinowitz_richter",
    "_stroud",
    "_stroud_secrest",
    #
    "schemes",
    "get_good_scheme",
]
