from . import _classical
from . import _ditkin
from . import _hammer_stroud
from . import _mysovskih
from . import _stroud

from ._helpers import schemes, get_good_scheme

__all__ = [
    "_classical",
    "_ditkin",
    "_hammer_stroud",
    "_mysovskih",
    "_stroud",
    #
    "schemes",
    "get_good_scheme",
]
